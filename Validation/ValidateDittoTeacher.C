///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   ValidateDittoTeacher.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Validate Ditto against its PYTHIA teacher.
///

// ValidateDittoTeacher.C
//
// Direct closure validation of Ditto against the PYTHIA teacher stored in
// the tune / generator card.
//
// Observables:
//   - P(Nch)
//   - <pT> vs Nch for all charged particles and every tuned species
//   - <N_species> vs Nch for every tuned species
//   - c2{2} vs Nch and v2{2} vs Nch from charged-particle Q-vectors
//
// Definitions:
//   Nch:
//     number of final charged particles in |eta| < tune.activityEtaMax,
//     exactly matching the tuner definition.
//
//   <pT>(Nch):
//     particle-weighted mean pT in |eta| < observableEtaMax.  For each Nch
//     class every accepted particle contributes one entry to the TProfile.
//
//   <N_species>(Nch):
//     event-averaged species yield in |eta| < observableEtaMax.  Each event
//     contributes once, including zero yield.
//
//   c2{2}(Nch):
//     pair-weighted two-particle azimuthal correlator
//
//       c2{2} = <cos[2(phi1 - phi2)]>
//
//     evaluated with charged particles in |eta| < observableEtaMax using
//
//       (|Q2|^2 - M) / [M(M - 1)].
//
//   v2{2} = sqrt(c2{2}) for positive c2{2}.  Negative bins are left empty.
//
// Important:
//   The current tune card declares azimuthModel = "uniform".  Therefore a
//   disagreement in c2{2}/v2{2} with PYTHIA is expected if PYTHIA contains
//   non-flow azimuthal correlations from jets, resonance decays, etc.  This
//   observable is deliberately included to expose that missing correlation.
//
// ROOT usage:
//   root -l
//   .L ValidateDittoTeacher.C+
//   ValidateDittoTeacher("Ditto_tune_pythia8_inel_136tev.root");
//
// The macro uses the complete PYTHIA card content persisted inside the tune,
// so the teacher sample is generated from the same physics configuration.

#include "Ditto.h"
#include "DittoTune.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace DittoTeacherValidation
{

struct SpeciesValidation {
  int pdg = 0;
  std::string name;

  TProfile* yieldTeacher = nullptr;
  TProfile* yieldFast = nullptr;

  TProfile* meanPtTeacher = nullptr;
  TProfile* meanPtFast = nullptr;
};

struct SampleHistograms {
  bool teacherSample = false;

  TH1D* nch = nullptr;

  TProfile* chargedMeanPt = nullptr;

  // Pair-weighted c2{2}; profile weight is M(M-1).
  TProfile* c2 = nullptr;

  std::vector<SpeciesValidation>* species = nullptr;
};

std::string trim(const std::string& input)
{
  const auto first = input.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return {};

  const auto last = input.find_last_not_of(" \t\r\n");
  return input.substr(first, last - first + 1);
}

void configureTeacherFromTune(Pythia8::Pythia& pythia,
                              const Ditto::Tune& tune,
                              int seed)
{
  std::istringstream card(tune.pythiaCardContent);
  std::string line;

  while (std::getline(card, line)) {
    const std::string stripped = trim(line);

    if (stripped.empty() ||
        stripped.front() == '!' ||
        stripped.front() == '#') {
      continue;
    }

    if (!pythia.readString(stripped)) {
      throw std::runtime_error("PYTHIA rejected tune-card line: " + stripped);
    }
  }

  // Generate an independent teacher validation sample rather than replaying
  // the random sequence used to build the tune.
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));

  if (!pythia.init()) {
    throw std::runtime_error("Could not initialize PYTHIA teacher from the tune card");
  }
}

std::vector<double> histogramEdges(const TH1& hist)
{
  std::vector<double> edges(hist.GetNbinsX() + 1);

  for (int i = 1; i <= hist.GetNbinsX(); ++i)
    edges[i - 1] = hist.GetXaxis()->GetBinLowEdge(i);

  edges.back() = hist.GetXaxis()->GetBinUpEdge(hist.GetNbinsX());

  return edges;
}

TProfile* makeProfile(const char* name,
                      const char* title,
                      const std::vector<double>& edges)
{
  auto* profile = new TProfile(name, title, static_cast<int>(edges.size()) - 1, edges.data());

  profile->SetDirectory(nullptr);
  return profile;
}

TH1D* profileToHistogram(const TProfile& profile,
                         const char* name,
                         const char* title)
{
  std::vector<double> edges(profile.GetNbinsX() + 1);

  for (int i = 1; i <= profile.GetNbinsX(); ++i)
    edges[i - 1] = profile.GetXaxis()->GetBinLowEdge(i);

  edges.back() = profile.GetXaxis()->GetBinUpEdge(profile.GetNbinsX());

  auto* hist = new TH1D(name, title, profile.GetNbinsX(), edges.data());
  hist->SetDirectory(nullptr);

  for (int i = 1; i <= profile.GetNbinsX(); ++i) {
    if (profile.GetBinEntries(i) <= 0.0)
      continue;

    hist->SetBinContent(i, profile.GetBinContent(i));
    hist->SetBinError(i, profile.GetBinError(i));
  }

  return hist;
}

TH1D* makeRatio(const TH1& fast,
                const TH1& teacher,
                const char* name,
                const char* title)
{
  auto* ratio = static_cast<TH1D*>(fast.Clone(name));

  ratio->Reset("ICES");
  ratio->SetDirectory(nullptr);
  ratio->SetTitle(title);

  for (int i = 1; i <= fast.GetNbinsX(); ++i) {
    const double f = fast.GetBinContent(i);
    const double fe = fast.GetBinError(i);
    const double t = teacher.GetBinContent(i);
    const double te = teacher.GetBinError(i);

    if (t == 0.0)
      continue;

    const double r = f / t;

    double error2 = 0.0;
    if (f != 0.0)
      error2 += (fe / f) * (fe / f);
    if (t != 0.0)
      error2 += (te / t) * (te / t);

    if ((r > 1.5 || r < 0.5) && 0) {
      std::cerr << "Warning: large ratio in bin " << i
                << ": fast = " << f
                << ", teacher = " << t
                << ", ratio = " << r
                << ", error = " << std::sqrt(error2) * r
                << std::endl;
    }
    ratio->SetBinContent(i, r);
    ratio->SetBinError(i, std::abs(r) * std::sqrt(error2));
  }

  return ratio;
}

TH1D* makeV2Histogram(const TProfile& c2,
                      const char* name,
                      const char* title)
{
  std::vector<double> edges(c2.GetNbinsX() + 1);

  for (int i = 1; i <= c2.GetNbinsX(); ++i)
    edges[i - 1] = c2.GetXaxis()->GetBinLowEdge(i);

  edges.back() = c2.GetXaxis()->GetBinUpEdge(c2.GetNbinsX());

  auto* v2 = new TH1D(name, title, c2.GetNbinsX(), edges.data());

  v2->SetDirectory(nullptr);

  for (int i = 1; i <= c2.GetNbinsX(); ++i) {
    if (c2.GetBinEntries(i) <= 0.0)
      continue;

    const double value = c2.GetBinContent(i);
    const double error = c2.GetBinError(i);

    if (value <= 0.0)
      continue;

    const double root = std::sqrt(value);
    v2->SetBinContent(i, root);
    v2->SetBinError(i, error / (2.0 * root));
  }

  return v2;
}

bool isFinalTeacher(const Pythia8::Particle& particle)
{
  return particle.isFinal();
}

bool isFinalFast(const Pythia8::Particle& particle)
{
  return particle.isFinal();
}

template <typename FinalSelector>
void analyzeEvent(const Pythia8::Event& event,
                  FinalSelector isFinal,
                  const Ditto::Tune& tune,
                  double observableEtaMax,
                  SampleHistograms& output)
{
  int nch = 0;
  int nNeutral = 0;
  int nNotInEta = 0;

  // First pass: exact activity multiplicity.
  for (int i = 1; i < event.size(); ++i) {
    const auto& particle = event[i];

    if (!isFinal(particle))
      continue;

    if (!particle.isCharged()) {
      nNeutral += 1;
      if (!particle.isNeutral()) {
        throw std::runtime_error("Not charged particle cannot be also not neutral");
      }
      continue;
    }

    if (std::abs(particle.eta()) > tune.activityEtaMax) {
      nNotInEta++;
      continue;
    }
    ++nch;
  }

  // if (nch <= 0) {
  //   std::cerr << "Warning: event has zero activity multiplicity (Nch = 0), neutral = " << nNeutral << ", outside eta range = " << nNotInEta << std::endl;
  // }
  output.nch->Fill(nch);

  std::vector<int> speciesCounts(output.species->size(), 0);

  std::complex<double> q2{0.0, 0.0};
  int qMultiplicity = 0;

  // Second pass: conditional observables.
  for (int i = 1; i < event.size(); ++i) {
    const auto& particle = event[i];

    if (!isFinal(particle))
      continue;

    const double eta = particle.eta();

    if (!std::isfinite(eta) || std::abs(eta) >= observableEtaMax) {
      continue;
    }

    const double pt = particle.pT();

    if (particle.isCharged()) {
      output.chargedMeanPt->Fill(nch, pt);

      const double phi = particle.phi();
      q2 += std::polar(1.0, 2.0 * phi);
      ++qMultiplicity;
    }

    const int pdg = particle.id();

    for (std::size_t iSpecies = 0; iSpecies < output.species->size(); ++iSpecies) {
      auto& species = (*output.species)[iSpecies];

      if (pdg != species.pdg)
        continue;

      ++speciesCounts[iSpecies];

      TProfile* meanPt = output.teacherSample
                           ? species.meanPtTeacher
                           : species.meanPtFast;

      meanPt->Fill(nch, pt);
      break;
    }
  }

  if (qMultiplicity >= 2) {
    const double pairs = static_cast<double>(qMultiplicity) * static_cast<double>(qMultiplicity - 1);

    const double numerator = std::norm(q2) - static_cast<double>(qMultiplicity);

    output.c2->Fill(nch, numerator / pairs, pairs);
  }

  // Fill event-wise species yields, including zeros.
  for (std::size_t iSpecies = 0; iSpecies < output.species->size(); ++iSpecies) {
    auto& species = (*output.species)[iSpecies];

    TProfile* yield = output.teacherSample
                        ? species.yieldTeacher
                        : species.yieldFast;

    yield->Fill(nch, speciesCounts[iSpecies]);
  }
}

void normalizeDistribution(TH1& hist)
{
  const double integral = hist.Integral();

  if (integral > 0.0)
    hist.Scale(1.0 / integral);
}

void styleTeacher(TH1& hist)
{
  hist.SetMarkerStyle(20);
  hist.SetMarkerSize(0.8);
  hist.SetLineWidth(2);
}

void styleFast(TH1& hist)
{
  hist.SetMarkerStyle(24);
  hist.SetMarkerSize(0.8);
  hist.SetLineWidth(2);
}

void drawWithRatio(TH1& teacher,
                   TH1& fast,
                   TH1& ratio,
                   const char* canvasName,
                   const char* legendDetail,
                   const std::string& pdfName,
                   bool first,
                   bool last,
                   bool logY = false)
{
  auto* canvas = new TCanvas(canvasName, canvasName, 800, 800);

  auto* upper = new TPad(Form("%s_upper", canvasName),
                         "",
                         0.0,
                         0.30,
                         1.0,
                         1.0);

  auto* lower = new TPad(Form("%s_lower", canvasName),
                         "",
                         0.0,
                         0.00,
                         1.0,
                         0.30);

  upper->SetBottomMargin(0.02);
  upper->SetLeftMargin(0.14);
  upper->SetRightMargin(0.04);
  upper->SetLogy(logY);

  lower->SetTopMargin(0.02);
  lower->SetBottomMargin(0.30);
  lower->SetLeftMargin(0.14);
  lower->SetRightMargin(0.04);

  upper->Draw();
  lower->Draw();

  upper->cd();

  styleTeacher(teacher);
  styleFast(fast);

  double maximum = std::max(teacher.GetMaximum(), fast.GetMaximum());

  if (maximum <= 0.0)
    maximum = 1.0;

  teacher.SetMaximum(logY ? maximum * 5.0 : maximum * 1.3);

  if (!logY)
    teacher.SetMinimum(0.0);

  teacher.GetXaxis()->SetLabelSize(0.0);
  teacher.GetXaxis()->SetTitleSize(0.0);

  teacher.Draw("E1");
  fast.Draw("E1 SAME");

  auto* legend = new TLegend(0.57, 0.70, 0.92, 0.89);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(&teacher, "PYTHIA teacher", "lep");
  legend->AddEntry(&fast, "Ditto", "lep");

  if (legendDetail && std::string(legendDetail).size() > 0)
    legend->AddEntry((TObject*)nullptr, legendDetail, "");

  legend->Draw();

  lower->cd();

  ratio.SetMarkerStyle(24);
  ratio.SetMarkerSize(0.7);
  ratio.GetYaxis()->SetRangeUser(0.5, 1.5);
  ratio.GetYaxis()->SetNdivisions(505);
  ratio.GetYaxis()->SetTitleSize(0.10);
  ratio.GetYaxis()->SetTitleOffset(0.55);
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

  if (first)
    printName += "(";
  else if (last)
    printName += ")";

  canvas->Print(printName.c_str());
}

void drawNoRatio(TH1& teacher,
                 TH1& fast,
                 const char* canvasName,
                 const char* legendDetail,
                 const std::string& pdfName,
                 bool first,
                 bool last)
{
  auto* canvas = new TCanvas(canvasName, canvasName, 800, 700);

  canvas->SetLeftMargin(0.14);
  canvas->SetRightMargin(0.04);

  styleTeacher(teacher);
  styleFast(fast);

  const double maximum = std::max(teacher.GetMaximum(), fast.GetMaximum());

  const double minimum = std::min(teacher.GetMinimum(), fast.GetMinimum());

  teacher.SetMaximum(maximum > 0.0 ? maximum * 1.3 : 1.0);

  if (minimum >= 0.0)
    teacher.SetMinimum(0.0);

  teacher.Draw("E1");
  fast.Draw("E1 SAME");

  auto* legend = new TLegend(0.57, 0.70, 0.92, 0.89);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(&teacher, "PYTHIA teacher", "lep");
  legend->AddEntry(&fast, "Ditto", "lep");

  if (legendDetail && std::string(legendDetail).size() > 0)
    legend->AddEntry((TObject*)nullptr, legendDetail, "");

  legend->Draw();

  std::string printName = pdfName;

  if (first)
    printName += "(";
  else if (last)
    printName += ")";

  canvas->Print(printName.c_str());
}

} // namespace DittoTeacherValidation

void ValidateDittoTeacher(const char* tuneFile = "Ditto_tune_pythia8_inel_136tev.root",
                          std::uint64_t nTeacherEvents = 1E5,
                          std::uint64_t nFastEvents = 1E5,
                          const char* outputFile = "Validation/Ditto_teacher_closure.root",
                          const char* outputPdf = "Validation/Ditto_teacher_closure.pdf",
                          int teacherSeed = 98765432,
                          std::uint64_t fastSeed = 123456789,
                          double observableEtaMax = -1.0)
{
  using namespace DittoTeacherValidation;

  gStyle->SetOptStat(0);

  auto tune = Ditto::Tune::load(tuneFile);

  if (!tune)
    throw std::runtime_error("Could not load Ditto tune");

  tune->validate(true);

  if (observableEtaMax <= 0.0)
    observableEtaMax = tune->activityEtaMax;

  if (observableEtaMax > tune->particleEtaMax) {
    throw std::invalid_argument("observableEtaMax cannot exceed the tune particleEtaMax");
  }

  const auto nchEdges = histogramEdges(tune->hNch);

  std::vector<SpeciesValidation> species;
  species.reserve(static_cast<std::size_t>(tune->numberOfSpecies()));

  for (int iSpecies = 0; iSpecies < tune->numberOfSpecies(); ++iSpecies) {
    const auto* entry = tune->speciesAt(iSpecies);

    if (!entry)
      continue;

    SpeciesValidation observable;
    observable.pdg = entry->pdg;
    observable.name = entry->particleName.empty() ? std::to_string(entry->pdg) : entry->particleName;

    observable.yieldTeacher = makeProfile(Form("pYieldVsNch_Teacher_PDG%d", entry->pdg),
                                          Form(";%s;<N_{%s}>", "N_{ch}", observable.name.c_str()),
                                          nchEdges);

    observable.yieldFast = makeProfile(Form("pYieldVsNch_Fast_PDG%d", entry->pdg),
                                       Form(";%s;<N_{%s}>", "N_{ch}", observable.name.c_str()),
                                       nchEdges);

    observable.meanPtTeacher = makeProfile(Form("pMeanPtVsNch_Teacher_PDG%d", entry->pdg),
                                           Form(";%s;<p_{T}>_{%s} (GeV/c)", "N_{ch}", observable.name.c_str()),
                                           nchEdges);

    observable.meanPtFast = makeProfile(Form("pMeanPtVsNch_Fast_PDG%d", entry->pdg),
                                        Form(";%s;<p_{T}>_{%s} (GeV/c)", "N_{ch}", observable.name.c_str()),
                                        nchEdges);

    species.push_back(std::move(observable));
  }

  auto* hNchTeacher = new TH1D("hNch_Teacher",
                               ";N_{ch};probability",
                               static_cast<int>(nchEdges.size()) - 1,
                               nchEdges.data());

  auto* hNchFast = new TH1D("hNch_Fast",
                            ";N_{ch};probability",
                            static_cast<int>(nchEdges.size()) - 1,
                            nchEdges.data());

  hNchTeacher->Sumw2();
  hNchFast->Sumw2();
  hNchTeacher->SetDirectory(nullptr);
  hNchFast->SetDirectory(nullptr);

  auto* pChargedMeanPtTeacher = makeProfile("pChargedMeanPtVsNch_Teacher",
                                            ";N_{ch};<#it{p}_{T}>_{ch} (GeV/c)",
                                            nchEdges);

  auto* pChargedMeanPtFast = makeProfile("pChargedMeanPtVsNch_Fast",
                                         ";N_{ch};<#it{p}_{T}>_{ch} (GeV/c)",
                                         nchEdges);

  auto* pC2Teacher = makeProfile("pC2VsNch_Teacher",
                                 ";N_{ch};c_{2}{2}",
                                 nchEdges);

  auto* pC2Fast = makeProfile("pC2VsNch_Fast",
                              ";N_{ch};c_{2}{2}",
                              nchEdges);

  SampleHistograms teacherOutput{true,
                                 hNchTeacher,
                                 pChargedMeanPtTeacher,
                                 pC2Teacher,
                                 &species};

  SampleHistograms fastOutput{false,
                              hNchFast,
                              pChargedMeanPtFast,
                              pC2Fast,
                              &species};

  // --------------------------------------------------------------------------
  // PYTHIA teacher.
  // --------------------------------------------------------------------------

  Pythia8::Pythia teacher;
  configureTeacherFromTune(teacher, *tune, teacherSeed);

  std::uint64_t teacherGenerated = 0;
  std::uint64_t teacherAttempts = 0;

  const std::uint64_t maxTeacherAttempts = std::max<std::uint64_t>(nTeacherEvents * 10, nTeacherEvents + 1000);

  while (teacherGenerated < nTeacherEvents) {
    if (++teacherAttempts > maxTeacherAttempts) {
      throw std::runtime_error("Too many failed PYTHIA events during validation");
    }

    if (!teacher.next())
      continue;

    analyzeEvent(teacher.event,
                 isFinalTeacher,
                 *tune,
                 observableEtaMax,
                 teacherOutput);

    ++teacherGenerated;

    if (teacherGenerated % 10000 == 0) {
      std::cout << "Teacher "
                << teacherGenerated
                << " / "
                << nTeacherEvents
                << "\n";
    }
  }

  // --------------------------------------------------------------------------
  // Ditto.
  // --------------------------------------------------------------------------

  Ditto::Config config;
  config.tuneFile = tuneFile;
  config.seed = fastSeed;
  config.enableTimingMetrics = false;
  config.enableDetailedTimingMetrics = false;

  Ditto::Generator generator(config);

  // The Event object used by loadParticles only needs a valid PYTHIA particle
  // data table; it does not generate events.
  Pythia8::Pythia fastEventOwner;
  fastEventOwner.event.init("Ditto event", &fastEventOwner.particleData);

  for (std::uint64_t iev = 0; iev < nFastEvents; ++iev) {
    generator.generate();
    generator.loadParticles(fastEventOwner.event, true);

    analyzeEvent(fastEventOwner.event,
                 isFinalFast,
                 *tune,
                 observableEtaMax,
                 fastOutput);

    if ((iev + 1) % 100000 == 0) {
      std::cout << "Ditto "
                << iev + 1
                << " / "
                << nFastEvents
                << "\n";
    }
  }

  // --------------------------------------------------------------------------
  // Derived histograms.
  // --------------------------------------------------------------------------

  normalizeDistribution(*hNchTeacher);
  normalizeDistribution(*hNchFast);

  auto* hRatioNch = makeRatio(*hNchFast,
                              *hNchTeacher,
                              "hRatio_Nch",
                              ";N_{ch};Ditto / PYTHIA");

  auto* hChargedMeanPtTeacher = profileToHistogram(*pChargedMeanPtTeacher,
                                                   "hChargedMeanPtVsNch_Teacher",
                                                   ";N_{ch};<#it{p}_{T}>_{ch} (GeV/c)");

  auto* hChargedMeanPtFast = profileToHistogram(*pChargedMeanPtFast,
                                                "hChargedMeanPtVsNch_Fast",
                                                ";N_{ch};<#it{p}_{T}>_{ch} (GeV/c)");

  auto* hRatioChargedMeanPt = makeRatio(*hChargedMeanPtFast,
                                        *hChargedMeanPtTeacher,
                                        "hRatio_ChargedMeanPt",
                                        ";N_{ch};Ditto / PYTHIA");

  auto* hC2Teacher = profileToHistogram(*pC2Teacher,
                                        "hC2VsNch_Teacher",
                                        ";N_{ch};c_{2}{2}");

  auto* hC2Fast = profileToHistogram(*pC2Fast,
                                     "hC2VsNch_Fast",
                                     ";N_{ch};c_{2}{2}");

  auto* hV2Teacher = makeV2Histogram(*pC2Teacher,
                                     "hV2VsNch_Teacher",
                                     ";N_{ch};v_{2}{2}");

  auto* hV2Fast = makeV2Histogram(*pC2Fast,
                                  "hV2VsNch_Fast",
                                  ";N_{ch};v_{2}{2}");

  // --------------------------------------------------------------------------
  // Save.
  // --------------------------------------------------------------------------

  std::unique_ptr<TFile> output(
    TFile::Open(
      outputFile,
      "RECREATE"));

  if (!output || output->IsZombie()) {
    throw std::runtime_error(std::string("Cannot create ") + outputFile);
  }

  output->cd();

  tune->Write("GeneratorCard");

  hNchTeacher->Write();
  hNchFast->Write();
  hRatioNch->Write();

  pChargedMeanPtTeacher->Write();
  pChargedMeanPtFast->Write();
  hChargedMeanPtTeacher->Write();
  hChargedMeanPtFast->Write();
  hRatioChargedMeanPt->Write();

  pC2Teacher->Write();
  pC2Fast->Write();
  hC2Teacher->Write();
  hC2Fast->Write();
  hV2Teacher->Write();
  hV2Fast->Write();

  for (auto& observable : species) {
    auto* directory = output->mkdir(Form("PDG%d", observable.pdg));

    directory->cd();

    observable.yieldTeacher->Write();
    observable.yieldFast->Write();
    observable.meanPtTeacher->Write();
    observable.meanPtFast->Write();

    output->cd();
  }

  // --------------------------------------------------------------------------
  // PDF.
  // --------------------------------------------------------------------------

  bool firstPage = true;

  drawWithRatio(*hNchTeacher,
                *hNchFast,
                *hRatioNch,
                "cNch",
                "activity multiplicity",
                outputPdf,
                firstPage,
                false,
                true);

  firstPage = false;

  drawWithRatio(*hChargedMeanPtTeacher,
                *hChargedMeanPtFast,
                *hRatioChargedMeanPt,
                "cChargedMeanPt",
                Form("|#eta| < %.2f", observableEtaMax),
                outputPdf,
                false,
                false);

  drawNoRatio(*hC2Teacher,
              *hC2Fast,
              "cC2",
              Form("charged, |#eta| < %.2f", observableEtaMax),
              outputPdf,
              false,
              false);

  drawNoRatio(*hV2Teacher,
              *hV2Fast,
              "cV2",
              Form("charged, |#eta| < %.2f", observableEtaMax),
              outputPdf,
              false,
              species.empty());

  for (std::size_t iSpecies = 0;
       iSpecies < species.size();
       ++iSpecies) {
    auto& observable = species[iSpecies];

    auto* hYieldTeacher = profileToHistogram(*observable.yieldTeacher,
                                             Form("hYieldVsNch_Teacher_PDG%d", observable.pdg),
                                             Form(";N_{ch};<N_{%s}>", observable.name.c_str()));

    auto* hYieldFast = profileToHistogram(*observable.yieldFast,
                                          Form("hYieldVsNch_Fast_PDG%d", observable.pdg),
                                          Form(";N_{ch};<N_{%s}>", observable.name.c_str()));

    auto* hYieldRatio = makeRatio(*hYieldFast,
                                  *hYieldTeacher,
                                  Form("hRatio_YieldVsNch_PDG%d", observable.pdg),
                                  ";N_{ch};Ditto / PYTHIA");

    auto* hMeanPtTeacher = profileToHistogram(*observable.meanPtTeacher,
                                              Form("hMeanPtVsNch_Teacher_PDG%d", observable.pdg),
                                              Form(";N_{ch};<p_{T}>_{%s} (GeV/c)", observable.name.c_str()));

    auto* hMeanPtFast = profileToHistogram(*observable.meanPtFast,
                                           Form("hMeanPtVsNch_Fast_PDG%d", observable.pdg),
                                           Form(";N_{ch};<p_{T}>_{%s} (GeV/c)", observable.name.c_str()));

    auto* hMeanPtRatio = makeRatio(*hMeanPtFast,
                                   *hMeanPtTeacher,
                                   Form("hRatio_MeanPtVsNch_PDG%d", observable.pdg),
                                   ";N_{ch};Ditto / PYTHIA");

    output->cd();

    auto* directory = dynamic_cast<TDirectory*>(output->Get(Form("PDG%d", observable.pdg)));

    if (directory) {
      directory->cd();
      hYieldTeacher->Write();
      hYieldFast->Write();
      hYieldRatio->Write();
      hMeanPtTeacher->Write();
      hMeanPtFast->Write();
      hMeanPtRatio->Write();
    }

    output->cd();

    drawWithRatio(*hYieldTeacher,
                  *hYieldFast,
                  *hYieldRatio,
                  Form("cYield_PDG%d", observable.pdg),
                  observable.name.c_str(),
                  outputPdf,
                  false,
                  false);

    const bool lastPage = iSpecies + 1 == species.size();

    drawWithRatio(*hMeanPtTeacher,
                  *hMeanPtFast,
                  *hMeanPtRatio,
                  Form("cMeanPt_PDG%d", observable.pdg),
                  observable.name.c_str(),
                  outputPdf,
                  false,
                  lastPage);
  }

  output->Close();

  std::cout << "\nTeacher closure validation complete.\n"
            << "  tune              : " << tuneFile << "\n"
            << "  sqrt(s_NN)        : " << tune->sqrtSNN << " GeV\n"
            << "  activity |eta|    : " << tune->activityEtaMax << "\n"
            << "  observable |eta|  : " << observableEtaMax << "\n"
            << "  PYTHIA events     : " << teacherGenerated << "\n"
            << "  Ditto events  : " << nFastEvents << "\n"
            << "  azimuth model     : " << tune->azimuthModel << "\n"
            << "  ROOT output       : " << outputFile << "\n"
            << "  plots             : " << outputPdf << "\n";
}
