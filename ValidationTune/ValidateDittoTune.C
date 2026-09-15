///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   ValidateDittoTune.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Validate the persisted Ditto tune.
///

#include "DittoTune.h"

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace DittoTuneValidation
{

struct ClosureMetrics {
  double yieldRatio = 0.0;
  double ratioRms = 0.0;
  double statisticalRms = 0.0;
  double excessRms = 0.0;
  double jsDistance = 0.0;
  double chi2Ndf = 0.0;
  double maxAbsDeviation = 0.0;
  double score = -1.0;
  int binsUsed = 0;
};

ClosureMetrics scoreNch;
ClosureMetrics scoreNSelected;

struct SpeciesValidation {
  const Ditto::TuneSpecies* tune = nullptr;

  TH1D* teacherMultiplicity = nullptr;
  TH1D* fastMultiplicity = nullptr;
  TH1D* ratioMultiplicity = nullptr;
  ClosureMetrics scoreMultiplicity;

  TH1D* teacherPt = nullptr;
  TH1D* fastPt = nullptr;
  TH1D* ratioPt = nullptr;
  ClosureMetrics scorePt;

  TH1D* teacherEta = nullptr;
  TH1D* fastEta = nullptr;
  TH1D* ratioEta = nullptr;
  ClosureMetrics scoreEta;
};

bool isCharged(int pdgCode)
{
  auto* particle = TDatabasePDG::Instance()->GetParticle(pdgCode);

  return particle != nullptr && std::abs(particle->Charge()) > 0.0;
}

TH1D* cloneEmpty(const TH1D& source,
                 const std::string& name,
                 const std::string& title)
{
  auto* histogram = static_cast<TH1D*>(source.Clone(name.c_str()));

  histogram->Reset("ICES");
  histogram->SetTitle(title.c_str());
  histogram->SetDirectory(nullptr);
  if (histogram->GetSumw2N() == 0)
    histogram->Sumw2();

  return histogram;
}

TH1D* projectionY(const TH2D& source,
                  const std::string& name,
                  const std::string& title)
{
  auto* histogram = source.ProjectionY(name.c_str(),
                                       1,
                                       source.GetNbinsX(),
                                       "e");

  histogram->SetTitle(title.c_str());
  histogram->SetDirectory(nullptr);

  return histogram;
}

void normalizeProbability(TH1D& histogram)
{
  const double integral = histogram.Integral(1, histogram.GetNbinsX());

  if (integral <= 0.0) {
    return;
  }

  histogram.Scale(1.0 / integral);
}

void normalizePerEventDensity(TH1D& histogram,
                              double nEvents)
{
  if (nEvents <= 0.0) {
    throw std::runtime_error("Ditto tune validation: cannot normalize to zero events");
  }

  for (int i = 1; i <= histogram.GetNbinsX(); ++i) {
    const double width = histogram.GetXaxis()->GetBinWidth(i);

    const double norm = nEvents * width;

    histogram.SetBinContent(i, histogram.GetBinContent(i) / norm);

    histogram.SetBinError(i, histogram.GetBinError(i) / norm);
  }
}

TH1D* makeRatio(const TH1D& fast,
                const TH1D& teacher,
                const std::string& name)
{
  auto* ratio = static_cast<TH1D*>(fast.Clone(name.c_str()));

  ratio->Reset("ICES");
  ratio->SetDirectory(nullptr);

  for (int i = 1; i <= fast.GetNbinsX(); ++i) {
    const double valueFast = fast.GetBinContent(i);
    const double errorFast = fast.GetBinError(i);

    const double valueTeacher = teacher.GetBinContent(i);
    const double errorTeacher = teacher.GetBinError(i);

    if (valueTeacher <= 0.0) {
      continue;
    }

    const double value = valueFast / valueTeacher;

    double relativeVariance = 0.0;

    if (valueFast > 0.0) {
      relativeVariance += errorFast * errorFast /
                          (valueFast * valueFast);
    }

    relativeVariance += errorTeacher * errorTeacher /
                        (valueTeacher * valueTeacher);

    ratio->SetBinContent(i, value);
    ratio->SetBinError(i,
                       value * std::sqrt(relativeVariance));
  }

  return ratio;
}

void restrictNchRange(TH1D& histogram,
                      double minimum,
                      double maximum)
{
  for (int i = 1; i <= histogram.GetNbinsX(); ++i) {
    const double center = histogram.GetXaxis()->GetBinCenter(i);

    if (center < minimum || center >= maximum) {
      histogram.SetBinContent(i, 0.0);
      histogram.SetBinError(i, 0.0);
    }
  }
}

double numberOfConditionedTeacherEvents(const Ditto::Tune& tune)
{
  double events = 0.0;

  for (int i = 1; i <= tune.hNch.GetNbinsX(); ++i) {
    const double center = tune.hNch.GetXaxis()->GetBinCenter(i);

    if (center < tune.activityEdges.front() ||
        center >= tune.activityEdges.back()) {
      continue;
    }

    events += tune.hNch.GetBinContent(i);
  }

  return events;
}

void styleTeacher(TH1D& histogram)
{
  histogram.SetMarkerStyle(20);
  histogram.SetMarkerSize(0.8);
  histogram.SetLineWidth(2);
}

void styleFast(TH1D& histogram)
{
  histogram.SetMarkerStyle(24);
  histogram.SetMarkerSize(0.8);
  histogram.SetLineWidth(2);
}

void drawComparison(const std::string& canvasName,
                    const std::string& label,
                    TH1D& fast,
                    TH1D& teacher,
                    TH1D& ratio,
                    const ClosureMetrics& metrics,
                    const std::string& pdfName,
                    bool first,
                    bool last,
                    bool logY,
                    double ratioMinimum = 0.5,
                    double ratioMaximum = 1.5)
{
  auto* canvas = new TCanvas(canvasName.c_str(),
                             label.c_str(),
                             800,
                             800);

  auto* upper = new TPad((canvasName + "_upper").c_str(),
                         "upper",
                         0.0,
                         0.30,
                         1.0,
                         1.0);

  auto* lower = new TPad((canvasName + "_lower").c_str(),
                         "lower",
                         0.0,
                         0.00,
                         1.0,
                         0.30);

  upper->SetBottomMargin(0.02);
  upper->SetLeftMargin(0.14);
  upper->SetRightMargin(0.04);

  if (logY) {
    upper->SetLogy();
  }

  lower->SetTopMargin(0.02);
  lower->SetBottomMargin(0.30);
  lower->SetLeftMargin(0.14);
  lower->SetRightMargin(0.04);

  upper->Draw();
  lower->Draw();

  upper->cd();

  styleTeacher(teacher);
  styleFast(fast);

  double yMaximum = std::max(fast.GetMaximum(), teacher.GetMaximum());

  if (yMaximum <= 0.0) {
    yMaximum = 1.0;
  }

  if (logY) {
    double yMinimum = std::numeric_limits<double>::max();

    for (int i = 1; i <= fast.GetNbinsX(); ++i) {
      const double valueFast = fast.GetBinContent(i);

      const double valueTeacher = teacher.GetBinContent(i);

      if (valueFast > 0.0) {
        yMinimum = std::min(yMinimum,
                            valueFast);
      }

      if (valueTeacher > 0.0) {
        yMinimum = std::min(yMinimum,
                            valueTeacher);
      }
    }

    if (!std::isfinite(yMinimum) ||
        yMinimum == std::numeric_limits<double>::max()) {
      yMinimum = 1.0e-8;
    }

    teacher.SetMinimum(yMinimum * 0.4);
    teacher.SetMaximum(yMaximum * 3.0);
  } else {
    teacher.SetMinimum(0.0);
    teacher.SetMaximum(yMaximum * 1.25);
  }

  teacher.GetXaxis()->SetLabelSize(0.0);
  teacher.GetXaxis()->SetTitleSize(0.0);
  teacher.GetYaxis()->SetTitleOffset(1.35);

  teacher.Draw("E1");
  fast.Draw("E1 SAME");

  auto* legend = new TLegend(0.52, 0.66, 0.91, 0.89);

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  legend->AddEntry(&teacher, "PYTHIA teacher", "lep");

  legend->AddEntry(&fast, "Ditto", "lep");

  legend->AddEntry((TObject*)nullptr, label.c_str(), "");

  if (metrics.score >= 0.0) {
    legend->AddEntry((TObject*)nullptr,
                     Form("closure score = %.1f / 100", metrics.score),
                     "");
  } else {
    legend->AddEntry((TObject*)nullptr,
                     "closure score = N/A",
                     "");
  }

  legend->Draw();

  lower->cd();

  ratio.SetMarkerStyle(24);
  ratio.SetMarkerSize(0.8);

  ratio.GetYaxis()->SetRangeUser(ratioMinimum,
                                 ratioMaximum);

  ratio.GetYaxis()->SetNdivisions(505);
  ratio.GetYaxis()->SetTitle("Ditto / PYTHIA");
  ratio.GetYaxis()->SetTitleSize(0.09);
  ratio.GetYaxis()->SetTitleOffset(0.62);
  ratio.GetYaxis()->SetLabelSize(0.09);

  ratio.GetXaxis()->SetTitleSize(0.12);
  ratio.GetXaxis()->SetLabelSize(0.10);

  ratio.Draw("E1");

  auto* unity = new TLine(ratio.GetXaxis()->GetXmin(),
                          1.0,
                          ratio.GetXaxis()->GetXmax(),
                          1.0);

  unity->SetLineStyle(2);
  unity->Draw("SAME");

  canvas->cd();

  std::string printName = pdfName;

  if (first) {
    printName += "(";
  } else if (last) {
    printName += ")";
  }

  canvas->Print(printName.c_str());
  canvas->Modified();
  canvas->Update();
}

void writeComparison(TDirectory* directory,
                     TH1D& fast,
                     TH1D& teacher,
                     TH1D& ratio)
{
  directory->cd();

  fast.Write("Ditto");
  teacher.Write("PythiaTeacher");
  ratio.Write("Ratio_Ditto_Pythia");
}

double densityIntegral(const TH1D& histogram)
{
  double integral = 0.0;

  for (int i = 1; i <= histogram.GetNbinsX(); ++i) {
    integral += histogram.GetBinContent(i) * histogram.GetXaxis()->GetBinWidth(i);
  }

  return integral;
}

double effectiveEntries(const TH1D& histogram,
                        int bin)
{
  const double value = histogram.GetBinContent(bin);
  const double error = histogram.GetBinError(bin);

  if (value <= 0.0 || error <= 0.0) {
    return 0.0;
  }

  return value * value / (error * error);
}

ClosureMetrics closureMetrics(const TH1D& fast,
                              const TH1D& teacher,
                              double minTeacherEffectiveEntries = 100.0,
                              double scoreScale = 0.05)
{
  ClosureMetrics metrics;

  const double yieldTeacher = densityIntegral(teacher);
  const double yieldFast = densityIntegral(fast);

  if (yieldTeacher > 0.0) {
    metrics.yieldRatio = yieldFast / yieldTeacher;
  }

  double ratioSquaredSum = 0.0;
  double statisticalVarianceSum = 0.0;
  double chi2 = 0.0;

  for (int i = 1; i <= teacher.GetNbinsX(); ++i) {
    const double valueTeacher = teacher.GetBinContent(i);

    if (valueTeacher <= 0.0) {
      continue;
    }

    if (effectiveEntries(teacher, i) < minTeacherEffectiveEntries) {
      continue;
    }

    const double valueFast = fast.GetBinContent(i);
    const double errorFast = fast.GetBinError(i);
    const double errorTeacher = teacher.GetBinError(i);

    const double ratio = valueFast / valueTeacher;
    const double deviation = ratio - 1.0;

    ratioSquaredSum += deviation * deviation;
    metrics.maxAbsDeviation = std::max(metrics.maxAbsDeviation,
                                       std::abs(deviation));

    // Propagated uncertainty for R = F / T.
    const double ratioVariance = (errorFast / valueTeacher) *
                                   (errorFast / valueTeacher) +
                                 (valueFast * errorTeacher /
                                  (valueTeacher * valueTeacher)) *
                                   (valueFast * errorTeacher /
                                    (valueTeacher * valueTeacher));

    statisticalVarianceSum += ratioVariance;

    const double variance = errorFast * errorFast + errorTeacher * errorTeacher;

    if (variance > 0.0) {
      const double difference = valueFast - valueTeacher;
      chi2 += difference * difference / variance;
    }

    ++metrics.binsUsed;
  }

  if (metrics.binsUsed > 0) {
    const double inverseBins = 1.0 / static_cast<double>(metrics.binsUsed);

    metrics.ratioRms = std::sqrt(ratioSquaredSum * inverseBins);
    metrics.statisticalRms = std::sqrt(statisticalVarianceSum * inverseBins);

    const double excessVariance = metrics.ratioRms * metrics.ratioRms - metrics.statisticalRms * metrics.statisticalRms;

    metrics.excessRms = std::sqrt(std::max(0.0, excessVariance));

    metrics.chi2Ndf = chi2 * inverseBins;

    if (scoreScale > 0.0) {
      const double scaled = metrics.excessRms / scoreScale;

      metrics.score = 100.0 * std::exp(-0.5 * scaled * scaled);
    }
  }

  if (yieldTeacher > 0.0 && yieldFast > 0.0) {
    std::vector<double> probabilityTeacher(teacher.GetNbinsX(), 0.0);
    std::vector<double> probabilityFast(fast.GetNbinsX(), 0.0);

    double normalizationTeacher = 0.0;
    double normalizationFast = 0.0;

    for (int i = 1; i <= teacher.GetNbinsX(); ++i) {
      const double width = teacher.GetXaxis()->GetBinWidth(i);

      probabilityTeacher[i - 1] = std::max(0.0, teacher.GetBinContent(i) * width);
      probabilityFast[i - 1] = std::max(0.0, fast.GetBinContent(i) * width);

      normalizationTeacher += probabilityTeacher[i - 1];
      normalizationFast += probabilityFast[i - 1];
    }

    if (normalizationTeacher > 0.0 && normalizationFast > 0.0) {
      double js = 0.0;

      for (int i = 0; i < teacher.GetNbinsX(); ++i) {
        const double p = probabilityTeacher[i] / normalizationTeacher;
        const double q = probabilityFast[i] / normalizationFast;
        const double m = 0.5 * (p + q);

        if (p > 0.0 && m > 0.0) {
          js += 0.5 * p * std::log(p / m);
        }

        if (q > 0.0 && m > 0.0) {
          js += 0.5 * q * std::log(q / m);
        }
      }

      metrics.jsDistance = std::sqrt(std::max(0.0, js) / std::log(2.0));
    }
  }

  return metrics;
}

void printClosureMetrics(const std::string& observable,
                         const ClosureMetrics& metrics)
{
  std::cout << std::left
            << std::setw(30)
            << observable;

  if (metrics.score < 0.0) {
    std::cout << "  score = N/A"
              << "  bins = "
              << metrics.binsUsed
              << "\n";
    return;
  }

  std::cout << std::right
            << std::fixed
            << std::setprecision(2)
            << "  score = "
            << std::setw(6)
            << metrics.score
            << " / 100"
            << std::setprecision(4)
            << "  ratio RMS = "
            << std::setw(8)
            << metrics.ratioRms
            << "  stat RMS = "
            << std::setw(8)
            << metrics.statisticalRms
            << "  excess RMS = "
            << std::setw(8)
            << metrics.excessRms
            << "  max |R-1| = "
            << std::setw(8)
            << metrics.maxAbsDeviation
            << "  JS = "
            << std::setw(8)
            << metrics.jsDistance
            << "  chi2/ndf = "
            << std::setw(8)
            << metrics.chi2Ndf
            << "  yield = "
            << std::setw(8)
            << metrics.yieldRatio
            << "  bins = "
            << metrics.binsUsed
            << "\n";
}

} // namespace DittoTuneValidation

void ValidateDittoTune(const char* inputFile = "Generation/Ditto.root",
                       const char* tuneFile = "Ditto_tune_pythia8_inel_136tev.root",
                      //  const char* tuneFile = "Ditto_tune_pythia8_PbPb_536tev.root",
                       const char* outputFile = "ValidationTune/Ditto_tune_validation.root",
                       const char* outputPdf = "ValidationTune/Ditto_vs_PythiaTune.pdf")
{
  using namespace DittoTuneValidation;

  gStyle->SetOptStat(0);

  auto tune = Ditto::Tune::load(tuneFile);

  const double nTeacherEvents = numberOfConditionedTeacherEvents(*tune);

  if (nTeacherEvents <= 0.0) {
    throw std::runtime_error("Ditto tune validation: no PYTHIA teacher events inside activity range");
  }

  std::unique_ptr<TFile> input(TFile::Open(inputFile, "READ"));

  if (!input || input->IsZombie()) {
    throw std::runtime_error(std::string("Cannot open Ditto file ") + inputFile);
  }

  auto* tree = dynamic_cast<TTree*>(input->Get("T"));

  if (!tree) {
    throw std::runtime_error("Could not find TTree 'T'");
  }

  TClonesArray* particles = nullptr;

  if (tree->SetBranchAddress("Particles", &particles) < 0) {
    throw std::runtime_error("Could not connect branch 'Particles'");
  }

  const Long64_t nEvents = tree->GetEntries();

  if (nEvents <= 0) {
    throw std::runtime_error("Ditto input contains no events");
  }

  std::cout << "Tune validation\n"
            << "  Ditto events       : "
            << nEvents
            << "\n"
            << "  PYTHIA teacher events  : "
            << nTeacherEvents
            << " inside activity range\n"
            << "  activity definition    : |eta| < "
            << tune->activityEtaMax
            << "\n"
            << "  particle acceptance    : |eta| < "
            << tune->particleEtaMax
            << "\n"
            << "  species                : "
            << tune->numberOfSpecies()
            << "\n";

  // -------------------------------------------------------------------------- // Event-level histograms.
  // --------------------------------------------------------------------------

  auto* teacherNch = static_cast<TH1D*>(tune->hNch.Clone("hTeacherNch"));

  teacherNch->SetDirectory(nullptr);

  restrictNchRange(*teacherNch, tune->activityEdges.front(), tune->activityEdges.back());

  auto* fastNch = cloneEmpty(*teacherNch, "hFastNch", ";N_{ch} (|#eta| < #eta_{activity});P(N_{ch})");

  teacherNch->SetTitle(";N_{ch} (|#eta| < #eta_{activity});P(N_{ch})");

  auto* teacherNSelected = static_cast<TH1D*>(tune->hNSelected.Clone("hTeacherNSelected"));

  teacherNSelected->SetDirectory(nullptr);
  teacherNSelected->SetTitle(";N_{selected};P(N_{selected})");

  auto* fastNSelected = cloneEmpty(*teacherNSelected, "hFastNSelected", ";N_{selected};P(N_{selected})");

  // -------------------------------------------------------------------------- // Species-level histograms.
  // --------------------------------------------------------------------------

  std::vector<SpeciesValidation> species;
  species.reserve(tune->numberOfSpecies());

  std::unordered_map<int, std::size_t> speciesIndex;
  speciesIndex.reserve(tune->numberOfSpecies());

  for (int i = 0; i < tune->numberOfSpecies(); ++i) {
    const auto* entry = tune->speciesAt(i);

    if (!entry) {
      continue;
    }

    SpeciesValidation validation;
    validation.tune = entry;

    const std::string tag = "pdg_" + std::to_string(entry->pdg);

    validation.teacherMultiplicity = projectionY(entry->hCountVsActivity,
                                                 "hTeacherMultiplicity_" + tag,
                                                 ";N_{species};P(N_{species})");

    validation.fastMultiplicity = cloneEmpty(*validation.teacherMultiplicity,
                                             "hFastMultiplicity_" + tag,
                                             ";N_{species};P(N_{species})");

    validation.teacherPt = projectionY(entry->hPtVsActivity,
                                       "hTeacherPt_" + tag,
                                       ";p_{T} (GeV/c);1/N_{ev} dN/dp_{T} ((GeV/c)^{-1})");

    validation.fastPt = cloneEmpty(*validation.teacherPt,
                                   "hFastPt_" + tag,
                                   ";p_{T} (GeV/c);1/N_{ev} dN/dp_{T} ((GeV/c)^{-1})");

    validation.teacherEta = projectionY(entry->hEtaVsActivity,
                                        "hTeacherEta_" + tag,
                                        ";#eta;1/N_{ev} dN/d#eta");

    validation.fastEta = cloneEmpty(*validation.teacherEta,
                                    "hFastEta_" + tag,
                                    ";#eta;1/N_{ev} dN/d#eta");

    speciesIndex[entry->pdg] = species.size();

    species.push_back(std::move(validation));
  }

  // -------------------------------------------------------------------------- // Ditto event loop.
  // --------------------------------------------------------------------------

  std::vector<int> speciesCounts(species.size(),
                                 0);

  Long64_t nchOverflowEvents = 0;
  Long64_t selectedOverflowEvents = 0;

  for (Long64_t iev = 0;
       iev < nEvents;
       ++iev) {
    tree->GetEntry(iev);

    if (!particles) {
      continue;
    }

    std::fill(speciesCounts.begin(),
              speciesCounts.end(),
              0);

    int nch = 0;
    int nSelected = 0;

    const int nParticles = particles->GetEntriesFast();

    for (int ip = 0;
         ip < nParticles;
         ++ip) {
      const auto* particle = static_cast<const TParticle*>(particles->UncheckedAt(ip));

      if (!particle) {
        continue;
      }

      if (particle->GetStatusCode() <= 0) {
        continue;
      }

      const int pdg = particle->GetPdgCode();

      const double eta = particle->Eta();

      if (std::isfinite(eta) &&
          std::abs(eta) < tune->activityEtaMax &&
          isCharged(pdg)) {
        ++nch;
      }

      const auto found = speciesIndex.find(pdg);

      if (found == speciesIndex.end()) {
        continue;
      }

      if (!std::isfinite(eta) ||
          std::abs(eta) >= tune->particleEtaMax) {
        continue;
      }

      const std::size_t iSpecies = found->second;

      ++speciesCounts[iSpecies];
      ++nSelected;

      species[iSpecies].fastPt->Fill(particle->Pt());

      species[iSpecies].fastEta->Fill(eta);
    }

    fastNch->Fill(nch);
    fastNSelected->Fill(nSelected);

    if (nch > fastNch->GetXaxis()->GetXmax()) {
      ++nchOverflowEvents;
    }

    if (nSelected >
        fastNSelected->GetXaxis()->GetXmax()) {
      ++selectedOverflowEvents;
    }

    for (std::size_t i = 0;
         i < species.size();
         ++i) {
      species[i].fastMultiplicity->Fill(speciesCounts[i]);
    }
  }

  // -------------------------------------------------------------------------- // Normalization and ratios.
  // --------------------------------------------------------------------------

  normalizeProbability(*teacherNch);
  normalizeProbability(*fastNch);

  normalizeProbability(*teacherNSelected);
  normalizeProbability(*fastNSelected);

  auto* ratioNch = makeRatio(*fastNch,
                             *teacherNch,
                             "hRatioNch");

  auto* ratioNSelected = makeRatio(*fastNSelected,
                                   *teacherNSelected,
                                   "hRatioNSelected");

  scoreNch = closureMetrics(*fastNch,
                            *teacherNch);

  scoreNSelected = closureMetrics(*fastNSelected,
                                  *teacherNSelected);

  for (auto& validation : species) {
    normalizeProbability(*validation.teacherMultiplicity);

    normalizeProbability(*validation.fastMultiplicity);

    normalizePerEventDensity(*validation.teacherPt,
                             nTeacherEvents);

    normalizePerEventDensity(*validation.fastPt,
                             static_cast<double>(nEvents));

    normalizePerEventDensity(*validation.teacherEta,
                             nTeacherEvents);

    normalizePerEventDensity(*validation.fastEta,
                             static_cast<double>(nEvents));

    const std::string tag = "pdg_" + std::to_string(validation.tune->pdg);

    validation.ratioMultiplicity = makeRatio(*validation.fastMultiplicity,
                                             *validation.teacherMultiplicity,
                                             "hRatioMultiplicity_" + tag);

    validation.ratioPt = makeRatio(*validation.fastPt,
                                   *validation.teacherPt,
                                   "hRatioPt_" + tag);

    validation.ratioEta = makeRatio(*validation.fastEta,
                                    *validation.teacherEta,
                                    "hRatioEta_" + tag);

    validation.scoreMultiplicity = closureMetrics(*validation.fastMultiplicity,
                                                  *validation.teacherMultiplicity);

    validation.scorePt = closureMetrics(*validation.fastPt,
                                        *validation.teacherPt);

    validation.scoreEta = closureMetrics(*validation.fastEta,
                                         *validation.teacherEta);
  }

  // -------------------------------------------------------------------------- // ROOT output.
  // --------------------------------------------------------------------------

  std::unique_ptr<TFile> output(TFile::Open(outputFile, "RECREATE"));

  if (!output || output->IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output file ") + outputFile);
  }

  auto* eventDirectory = output->mkdir("Event");

  auto* nchDirectory = eventDirectory->mkdir("Nch");

  writeComparison(nchDirectory,
                  *fastNch,
                  *teacherNch,
                  *ratioNch);

  auto* selectedDirectory = eventDirectory->mkdir("NSelected");

  writeComparison(selectedDirectory,
                  *fastNSelected,
                  *teacherNSelected,
                  *ratioNSelected);

  auto* speciesDirectory = output->mkdir("Species");

  for (auto& validation : species) {
    speciesDirectory->cd();

    const std::string directoryName = "pdg_" + std::to_string(validation.tune->pdg);

    auto* directory = speciesDirectory->mkdir(directoryName.c_str());

    auto* multiplicityDirectory = directory->mkdir("Multiplicity");

    writeComparison(multiplicityDirectory,
                    *validation.fastMultiplicity,
                    *validation.teacherMultiplicity,
                    *validation.ratioMultiplicity);

    auto* ptDirectory = directory->mkdir("Pt");

    writeComparison(ptDirectory,
                    *validation.fastPt,
                    *validation.teacherPt,
                    *validation.ratioPt);

    auto* etaDirectory = directory->mkdir("Eta");

    writeComparison(etaDirectory,
                    *validation.fastEta,
                    *validation.teacherEta,
                    *validation.ratioEta);
  }

  // -------------------------------------------------------------------------- // PDF.
  // --------------------------------------------------------------------------

  const int numberOfPages = 2 + 3 * static_cast<int>(species.size());

  int page = 0;

  drawComparison("c_Nch",
                 "event charged multiplicity",
                 *fastNch,
                 *teacherNch,
                 *ratioNch,
                 scoreNch,
                 outputPdf,
                 page == 0,
                 page == numberOfPages - 1,
                 true,
                 0.5,
                 1.5);

  ++page;

  drawComparison("c_NSelected",
                 "selected-particle multiplicity",
                 *fastNSelected,
                 *teacherNSelected,
                 *ratioNSelected,
                 scoreNSelected,
                 outputPdf,
                 page == 0,
                 page == numberOfPages - 1,
                 true,
                 0.5,
                 1.5);

  ++page;

  for (auto& validation : species) {
    const std::string particleLabel = validation.tune->particleName.empty()
                                        ? std::to_string(validation.tune->pdg)
                                        : validation.tune->particleName;

    drawComparison("c_Multiplicity_" + std::to_string(validation.tune->pdg),
                   particleLabel + " multiplicity",
                   *validation.fastMultiplicity,
                   *validation.teacherMultiplicity,
                   *validation.ratioMultiplicity,
                   validation.scoreMultiplicity,
                   outputPdf,
                   page == 0,
                   page == numberOfPages - 1,
                   true,
                   0.5,
                   1.5);

    ++page;

    drawComparison("c_Pt_" + std::to_string(validation.tune->pdg),
                   particleLabel + " p_{T}",
                   *validation.fastPt,
                   *validation.teacherPt,
                   *validation.ratioPt,
                   validation.scorePt,
                   outputPdf,
                   page == 0,
                   page == numberOfPages - 1,
                   true,
                   0.5,
                   1.5);

    ++page;

    drawComparison("c_Eta_" + std::to_string(validation.tune->pdg),
                   particleLabel + " #eta",
                   *validation.fastEta,
                   *validation.teacherEta,
                   *validation.ratioEta,
                   validation.scoreEta,
                   outputPdf,
                   page == 0,
                   page == numberOfPages - 1,
                   false,
                   0.5,
                   1.5);

    ++page;
  }

  output->cd();

  auto* scoreTree = new TTree("ClosureScores",
                              "Ditto tune closure scores");

  std::string scoreObservable;
  int scorePdg = 0;
  std::string scoreType;
  double scoreValue = 0.0;
  double scoreYieldRatio = 0.0;
  double scoreRatioRms = 0.0;
  double scoreStatisticalRms = 0.0;
  double scoreExcessRms = 0.0;
  double scoreJsDistance = 0.0;
  double scoreChi2Ndf = 0.0;
  double scoreMaxAbsDeviation = 0.0;
  int scoreBinsUsed = 0;

  scoreTree->Branch("observable", &scoreObservable);
  scoreTree->Branch("pdg", &scorePdg);
  scoreTree->Branch("type", &scoreType);
  scoreTree->Branch("score", &scoreValue);
  scoreTree->Branch("yieldRatio", &scoreYieldRatio);
  scoreTree->Branch("ratioRms", &scoreRatioRms);
  scoreTree->Branch("statisticalRms", &scoreStatisticalRms);
  scoreTree->Branch("excessRms", &scoreExcessRms);
  scoreTree->Branch("jsDistance", &scoreJsDistance);
  scoreTree->Branch("chi2Ndf", &scoreChi2Ndf);
  scoreTree->Branch("maxAbsDeviation", &scoreMaxAbsDeviation);
  scoreTree->Branch("binsUsed", &scoreBinsUsed);

  const auto fillScore = [&](const std::string& observable,
                             int pdg,
                             const std::string& type,
                             const ClosureMetrics& metrics) {
    scoreObservable = observable;
    scorePdg = pdg;
    scoreType = type;
    scoreValue = metrics.score;
    scoreYieldRatio = metrics.yieldRatio;
    scoreRatioRms = metrics.ratioRms;
    scoreStatisticalRms = metrics.statisticalRms;
    scoreExcessRms = metrics.excessRms;
    scoreJsDistance = metrics.jsDistance;
    scoreChi2Ndf = metrics.chi2Ndf;
    scoreMaxAbsDeviation = metrics.maxAbsDeviation;
    scoreBinsUsed = metrics.binsUsed;
    scoreTree->Fill();
  };

  fillScore("Nch", 0, "event", scoreNch);
  fillScore("NSelected", 0, "event", scoreNSelected);

  for (const auto& validation : species) {
    const std::string particleLabel = validation.tune->particleName.empty()
                                        ? std::to_string(validation.tune->pdg)
                                        : validation.tune->particleName;

    fillScore(particleLabel + " multiplicity",
              validation.tune->pdg,
              "multiplicity",
              validation.scoreMultiplicity);

    fillScore(particleLabel + " pT",
              validation.tune->pdg,
              "pt",
              validation.scorePt);

    fillScore(particleLabel + " eta",
              validation.tune->pdg,
              "eta",
              validation.scoreEta);
  }

  scoreTree->Write();

  output->Write();
  output->Close();

  // -------------------------------------------------------------------------- // Console summary.
  // --------------------------------------------------------------------------

  std::cout << "\nClosure scores\n"
            << "  ratio RMS uses the full Ditto/PYTHIA ratio.\n"
            << "  stat RMS is the expected RMS from propagated histogram statistics.\n"
            << "  excess RMS = sqrt(max(ratioRMS^2 - statRMS^2, 0)).\n"
            << "  Bins with fewer than 100 effective PYTHIA entries are excluded.\n"
            << "  Mapping: score = 100 exp[-0.5 (excessRMS / 0.05)^2].\n\n";

  printClosureMetrics("Nch", scoreNch);
  printClosureMetrics("NSelected", scoreNSelected);

  for (const auto& validation : species) {
    const std::string particleLabel = validation.tune->particleName.empty()
                                        ? std::to_string(validation.tune->pdg)
                                        : validation.tune->particleName;

    printClosureMetrics(particleLabel + " multiplicity", validation.scoreMultiplicity);

    printClosureMetrics(particleLabel + " pT", validation.scorePt);

    printClosureMetrics(particleLabel + " eta", validation.scoreEta);
  }

  std::cout << "\nMean / yield cross-checks\n"
            << "  <Nch> PYTHIA     : "
            << teacherNch->GetMean()
            << "\n"
            << "  <Nch> Ditto  : "
            << fastNch->GetMean()
            << "\n"
            << "  ratio            : "
            << (teacherNch->GetMean() > 0.0
                  ? fastNch->GetMean() /
                      teacherNch->GetMean()
                  : 0.0)
            << "\n"
            << "  <Nselected> PYTHIA    : "
            << teacherNSelected->GetMean()
            << "\n"
            << "  <Nselected> Ditto : "
            << fastNSelected->GetMean()
            << "\n";

  for (const auto& validation : species) {
    const double yieldTeacher = densityIntegral(*validation.teacherPt);
    const double yieldFast = densityIntegral(*validation.fastPt);

    std::cout << "  PDG "
              << std::setw(6)
              << validation.tune->pdg
              << "  <N> PYTHIA = "
              << std::setw(10)
              << yieldTeacher
              << "  Ditto = "
              << std::setw(10)
              << yieldFast
              << "  ratio = "
              << (yieldTeacher > 0.0
                    ? yieldFast / yieldTeacher
                    : 0.0)
              << "\n";
  }

  if (nchOverflowEvents > 0) {
    std::cout << "  WARNING Ditto Nch overflow events: "
              << nchOverflowEvents
              << "\n";
  }

  if (selectedOverflowEvents > 0) {
    std::cout << "  WARNING Ditto Nselected overflow events: "
              << selectedOverflowEvents
              << "\n";
  }

  std::cout << "\nTune validation complete.\n"
            << "  ROOT output : "
            << outputFile
            << "\n"
            << "  plots       : "
            << outputPdf
            << "\n";
}
