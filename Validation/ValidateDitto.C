///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   ValidateDitto.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Validate Ditto against identified-hadron spectra.
///

// ValidateDitto.C
//
// Compare Ditto TTree output with ALICE identified-hadron spectra in
// inelastic pp collisions at sqrt(s) = 13 TeV, |y| < 0.5.
//
// HEPData record: ins1797443 / DOI 10.17182/hepdata.100303
//   Table 1 : pi+ + pi-
//   Table 2 : K+  + K-
//   Table 6 : p   + pbar
//
// Charged-particle pseudorapidity reference:
//   HEPData record: ins1395253 / DOI 10.17182/hepdata.70847
//   Table 1, dependent variable y1: INEL dN_ch/deta at sqrt(s)=13 TeV
//
// The published identified-hadron observable is
//
//   (1/N_INEL) d^2N / (dpT dy)   [(GeV/c)^-1]
//
// so the MC bin content is normalized by
//
//   N_events * Delta y * Delta pT
//
// with Delta y = 1 for |y| < 0.5.
//
// ROOT usage:
//   root -l
//   .L ValidateDitto.C+
//   ValidateDitto("Ditto.root");
//
// On first use, the macro downloads the three HEPData ROOT tables with curl
// into the directory given by hepDataDir.  Set downloadHEPData=false if the
// files are already present or if the machine has no network access.

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TParticle.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace DittoValidation
{

struct ReferenceSpec {
  const char* shortName;
  const char* label;
  int absPdg;
  const char* tableName;
  int color;
  int marker;
};

static const ReferenceSpec kReferences[] = {
  {"Pion", "#pi^{+} + #pi^{-}", 211, "Table1", kRed + 1, 20},
  {"Kaon", "K^{+} + K^{-}", 321, "Table2", kGreen + 2, 21},
  {"Proton", "p + #bar{p}", 2212, "Table6", kBlue + 1, 22},
};

constexpr double kYMax = 0.5;
constexpr double kDeltaY = 2.0 * kYMax;
constexpr const char* kHepDataRecord = "ins1797443";
constexpr int kHepDataVersion = 1;

// ALICE charged-particle pseudorapidity distribution at sqrt(s)=13 TeV.
// HEPData 70847, Table 1 has two dependent variables:
//   y1 = INEL, y2 = INEL>0.
// Ditto currently normalizes to every generated event, so we compare to INEL.
constexpr const char* kChargedEtaHepDataRecord = "ins1395253";
constexpr int kChargedEtaHepDataVersion = 1;
constexpr const char* kChargedEtaTableName = "Table1";
constexpr int kChargedEtaYIndex = 1;

std::string hepDataFileName(const std::string& dir,
                            const ReferenceSpec& ref)
{
  return dir + "/HEPData-" + kHepDataRecord + "-v" +
         std::to_string(kHepDataVersion) + "-" + ref.tableName + ".root";
}

bool fileExists(const std::string& path)
{
  return gSystem->AccessPathName(path.c_str()) == kFALSE;
}

void downloadTable(const std::string& dir,
                   const ReferenceSpec& ref)
{
  gSystem->mkdir(dir.c_str(), true);

  const std::string output = hepDataFileName(dir, ref);
  if (fileExists(output))
    return;

  // HEPData documented record export endpoint.  Table1 without a space is
  // accepted by HEPData and avoids shell escaping issues.
  const std::string url = std::string("https://www.hepdata.net/record/") + kHepDataRecord + "?format=root&table=" + ref.tableName + "&version=" + std::to_string(kHepDataVersion);

  const std::string command = "curl -L --fail --silent --show-error '" + url + "' -o '" + output + "'";

  std::cout << "Downloading HEPData " << ref.tableName << " ...\n";
  const int rc = gSystem->Exec(command.c_str());

  if (rc != 0 || !fileExists(output)) {
    throw std::runtime_error("Could not download " + url + ". Download the table manually or rerun with network access.");
  }
}

std::string chargedEtaHepDataFileName(const std::string& dir)
{
  return dir + "/HEPData-" + kChargedEtaHepDataRecord + "-v" +
         std::to_string(kChargedEtaHepDataVersion) + "-" +
         kChargedEtaTableName + ".root";
}

void downloadChargedEtaTable(const std::string& dir)
{
  gSystem->mkdir(dir.c_str(), true);

  const std::string output = chargedEtaHepDataFileName(dir);
  if (fileExists(output))
    return;

  const std::string url =
    std::string("https://www.hepdata.net/record/") +
    kChargedEtaHepDataRecord + "?format=root&table=" +
    kChargedEtaTableName + "&version=" +
    std::to_string(kChargedEtaHepDataVersion);

  const std::string command =
    "curl -L --fail --silent --show-error '" + url +
    "' -o '" + output + "'";

  std::cout << "Downloading HEPData charged-particle dN/deta "
            << kChargedEtaTableName << " ...\n";
  const int rc = gSystem->Exec(command.c_str());

  if (rc != 0 || !fileExists(output)) {
    throw std::runtime_error(
      "Could not download " + url +
      ". Download the table manually or rerun with network access.");
  }
}

template <typename T>
T* findFirstObject(TDirectory* dir)
{
  if (!dir)
    return nullptr;

  TIter next(dir->GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    TObject* obj = key->ReadObj();
    if (!obj)
      continue;

    if (auto* wanted = dynamic_cast<T*>(obj))
      return wanted;

    if (auto* subdir = dynamic_cast<TDirectory*>(obj)) {
      if (auto* wanted = findFirstObject<T>(subdir))
        return wanted;
    }
  }

  return nullptr;
}

struct HepDataObjects {
  std::unique_ptr<TFile> file;
  TH1* hist = nullptr;
  TGraphAsymmErrors* graph = nullptr;
};

HepDataObjects openHEPData(const std::string& fileName,
                           const ReferenceSpec& ref)
{
  HepDataObjects out;
  out.file.reset(TFile::Open(fileName.c_str(), "READ"));

  if (!out.file || out.file->IsZombie())
    throw std::runtime_error("Cannot open HEPData file " + fileName);

  // The HEPData ROOT converter normally stores a single table in a directory
  // called "Table 1", "Table 2", ... and the central histogram / graph as
  // Hist1D_y1 / Graph1D_y1.  Try the canonical path first and fall back to a
  // recursive search to remain robust against converter naming changes.
  std::string spaced = ref.tableName;
  const auto pos = spaced.find("Table");
  if (pos == 0)
    spaced.insert(5, " ");

  out.hist = dynamic_cast<TH1*>(out.file->Get((spaced + "/Hist1D_y1").c_str()));
  out.graph = dynamic_cast<TGraphAsymmErrors*>(out.file->Get((spaced + "/Graph1D_y1").c_str()));

  if (!out.hist)
    out.hist = findFirstObject<TH1>(out.file.get());

  if (!out.graph)
    out.graph = findFirstObject<TGraphAsymmErrors>(out.file.get());

  if (!out.hist) {
    throw std::runtime_error("Could not find a 1D HEPData histogram in " + fileName);
  }

  out.hist->SetDirectory(nullptr);

  return out;
}

HepDataObjects openChargedEtaHEPData(const std::string& fileName)
{
  HepDataObjects out;
  out.file.reset(TFile::Open(fileName.c_str(), "READ"));

  if (!out.file || out.file->IsZombie())
    throw std::runtime_error("Cannot open HEPData file " + fileName);

  const std::string dirName = "Table 1";
  const std::string histName =
    dirName + "/Hist1D_y" + std::to_string(kChargedEtaYIndex);
  const std::string graphName =
    dirName + "/Graph1D_y" + std::to_string(kChargedEtaYIndex);

  out.hist = dynamic_cast<TH1*>(out.file->Get(histName.c_str()));
  out.graph = dynamic_cast<TGraphAsymmErrors*>(out.file->Get(graphName.c_str()));

  if (!out.hist) {
    throw std::runtime_error(
      "Could not find INEL charged-particle HEPData histogram '" +
      histName + "' in " + fileName);
  }

  // Detach the histogram from the input file so it stays valid while the
  // HepDataObjects object owns/closes the file.  The graph remains file-owned
  // for the lifetime of HepDataObjects.
  out.hist->SetDirectory(nullptr);
  return out;
}

std::vector<double> binEdges(const TH1& h)
{
  const int nBins = h.GetNbinsX();
  std::vector<double> edges(nBins + 1);

  for (int i = 1; i <= nBins; ++i)
    edges[i - 1] = h.GetXaxis()->GetBinLowEdge(i);

  edges[nBins] = h.GetXaxis()->GetBinUpEdge(nBins);
  return edges;
}

double rapidity(const TParticle& p)
{
  const double ep = p.Energy() + p.Pz();
  const double em = p.Energy() - p.Pz();

  if (ep <= 0.0 || em <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();

  return 0.5 * std::log(ep / em);
}

TH1D* makeMCHistogram(const TH1& data,
                      const ReferenceSpec& ref)
{
  const auto edges = binEdges(data);
  auto h = new TH1D(Form("hMC_%s", ref.shortName), Form(";%s;%s", "p_{T} (GeV/c)", "1/N_{INEL} d^{2}N/(dp_{T}dy) ((GeV/c)^{-1})"), static_cast<int>(edges.size()) - 1, edges.data());

  h->Sumw2();
  h->SetDirectory(nullptr);
  return h;
}

TH1D* makeMCEtaHistogram(const TH1& data)
{
  const auto edges = binEdges(data);
  auto* h = new TH1D(
    "hMC_Charged_Eta",
    ";#eta;1/N_{INEL} dN_{ch}/d#eta",
    static_cast<int>(edges.size()) - 1,
    edges.data());

  h->Sumw2();
  h->SetDirectory(nullptr);
  return h;
}

void normalizeSpectrum(TH1& h, Long64_t nEvents)
{
  if (nEvents <= 0)
    throw std::runtime_error("Cannot normalize to zero events");

  for (int i = 1; i <= h.GetNbinsX(); ++i) {
    const double width = h.GetXaxis()->GetBinWidth(i);
    const double norm = static_cast<double>(nEvents) * kDeltaY * width;

    h.SetBinContent(i, h.GetBinContent(i) / norm);
    h.SetBinError(i, h.GetBinError(i) / norm);
  }
}

TH1D* makeRatio(const TH1D& mc,
                const TH1& data,
                const ReferenceSpec& ref)
{
  auto ratio = static_cast<TH1D*>(mc.Clone(Form("hRatio_%s", ref.shortName)));
  ratio->Reset("ICES");
  ratio->SetDirectory(nullptr);
  ratio->SetTitle(";p_{T} (GeV/c);MC / data");

  for (int i = 1; i <= mc.GetNbinsX(); ++i) {
    const double m = mc.GetBinContent(i);
    const double me = mc.GetBinError(i);
    const double d = data.GetBinContent(i);

    if (d <= 0.0)
      continue;

    ratio->SetBinContent(i, m / d);

    // MC statistical uncertainty only. Experimental errors are shown on the
    // HEPData points in the upper panel and are intentionally not folded into
    // this histogram error.
    ratio->SetBinError(i, me / d);
  }

  return ratio;
}

TH1D* makeEtaRatio(const TH1D& mc, const TH1& data)
{
  auto* ratio = static_cast<TH1D*>(mc.Clone("hRatio_Charged_Eta"));
  ratio->Reset("ICES");
  ratio->SetDirectory(nullptr);
  ratio->SetTitle(";#eta;MC / data");

  for (int i = 1; i <= mc.GetNbinsX(); ++i) {
    const double m = mc.GetBinContent(i);
    const double me = mc.GetBinError(i);
    const double d = data.GetBinContent(i);

    if (d <= 0.0)
      continue;

    ratio->SetBinContent(i, m / d);
    ratio->SetBinError(i, me / d);
  }

  return ratio;
}

void styleReference(TGraphAsymmErrors* graph,
                    const ReferenceSpec& ref)
{
  if (!graph)
    return;

  graph->SetMarkerStyle(ref.marker);
  graph->SetMarkerSize(1.0);
  graph->SetMarkerColor(ref.color);
  graph->SetLineColor(ref.color);
}

void styleMC(TH1& h)
{
  h.SetMarkerStyle(24);
  h.SetMarkerSize(0.9);
  h.SetLineWidth(2);
}

void drawComparison(const ReferenceSpec& ref,
                    TH1D& mc,
                    TH1& dataHist,
                    TGraphAsymmErrors* dataGraph,
                    TH1D& ratio,
                    const std::string& pdfName,
                    bool first,
                    bool last)
{
  auto* canvas = new TCanvas(Form("c_%s", ref.shortName), ref.label, 800, 800);

  auto* upper = new TPad("upper", "upper", 0.0, 0.30, 1.0, 1.0);
  auto* lower = new TPad("lower", "lower", 0.0, 0.00, 1.0, 0.30);

  upper->SetBottomMargin(0.02);
  upper->SetLeftMargin(0.14);
  upper->SetRightMargin(0.04);
  upper->SetLogy();

  lower->SetTopMargin(0.02);
  lower->SetBottomMargin(0.30);
  lower->SetLeftMargin(0.14);
  lower->SetRightMargin(0.04);

  upper->Draw();
  lower->Draw();

  upper->cd();

  styleMC(mc);
  styleReference(dataGraph, ref);

  double yMin = 1e30;
  double yMax = 0.0;
  for (int i = 1; i <= mc.GetNbinsX(); ++i) {
    const double m = mc.GetBinContent(i);
    const double d = dataHist.GetBinContent(i);
    if (m > 0.0)
      yMin = std::min(yMin, m);
    if (d > 0.0)
      yMin = std::min(yMin, d);
    yMax = std::max(yMax, std::max(m, d));
  }

  if (!std::isfinite(yMin) || yMin <= 0.0)
    yMin = 1e-8;
  if (yMax <= 0.0)
    yMax = 1.0;

  mc.SetMinimum(yMin * 0.4);
  mc.SetMaximum(yMax * 3.0);
  mc.GetXaxis()->SetLabelSize(0.0);
  mc.GetXaxis()->SetTitleSize(0.0);
  mc.GetYaxis()->SetTitleOffset(1.35);
  mc.Draw("E1");

  if (dataGraph)
    dataGraph->Draw("P SAME");
  else {
    dataHist.SetMarkerStyle(ref.marker);
    dataHist.SetMarkerColor(ref.color);
    dataHist.SetLineColor(ref.color);
    dataHist.Draw("E1 SAME");
  }

  auto* legend = new TLegend(0.56, 0.69, 0.92, 0.89);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(&mc, "Ditto", "lep");
  if (dataGraph)
    legend->AddEntry(dataGraph, "ALICE, pp #sqrt{s}=13 TeV", "lep");
  else
    legend->AddEntry(&dataHist, "ALICE, pp #sqrt{s}=13 TeV", "lep");
  legend->AddEntry((TObject*)nullptr, ref.label, "");
  legend->AddEntry((TObject*)nullptr, "|y| < 0.5, INEL", "");
  legend->Draw();

  lower->cd();
  ratio.SetMarkerStyle(24);
  ratio.SetMarkerSize(0.8);
  ratio.GetYaxis()->SetRangeUser(0.0, 2.0);
  ratio.GetYaxis()->SetNdivisions(505);
  ratio.GetYaxis()->SetTitleSize(0.10);
  ratio.GetYaxis()->SetTitleOffset(0.55);
  ratio.GetYaxis()->SetLabelSize(0.09);
  ratio.GetXaxis()->SetTitleSize(0.12);
  ratio.GetXaxis()->SetLabelSize(0.10);
  ratio.Draw("E1");

  auto* unity = new TLine(ratio.GetXaxis()->GetXmin(), 1.0,
                          ratio.GetXaxis()->GetXmax(), 1.0);
  unity->SetLineStyle(2);
  unity->Draw("SAME");

  canvas->cd();
  canvas->Modified();
  canvas->Update();
  std::string printName = pdfName;
  if (first)
    printName += "(";
  else if (last)
    printName += ")";
  canvas->Print(printName.c_str());
}

bool isCharged(int pdgCode)
{
  auto* pdg = TDatabasePDG::Instance()->GetParticle(pdgCode);
  return (pdg != nullptr && std::abs(pdg->Charge()) > 0.0);
}

void normalizeEtaDistribution(TH1& h, Long64_t nEvents)
{
  if (nEvents <= 0)
    throw std::runtime_error("Cannot normalize to zero events");

  for (int i = 1; i <= h.GetNbinsX(); ++i) {
    const double width = h.GetXaxis()->GetBinWidth(i);
    const double norm = static_cast<double>(nEvents) * width;

    h.SetBinContent(i, h.GetBinContent(i) / norm);
    h.SetBinError(i, h.GetBinError(i) / norm);
  }
}

void drawChargedEta(TH1D& mc,
                    TH1& dataHist,
                    TGraphAsymmErrors* dataGraph,
                    TH1D& ratio,
                    const std::string& pdfName,
                    bool first,
                    bool last)
{
  auto* canvas = new TCanvas("c_Charged_Eta", "Charged Particles vs Eta", 800, 800);

  auto* upper = new TPad("upperEta", "upperEta", 0.0, 0.30, 1.0, 1.0);
  auto* lower = new TPad("lowerEta", "lowerEta", 0.0, 0.00, 1.0, 0.30);

  upper->SetBottomMargin(0.02);
  upper->SetLeftMargin(0.14);
  upper->SetRightMargin(0.04);

  lower->SetTopMargin(0.02);
  lower->SetBottomMargin(0.30);
  lower->SetLeftMargin(0.14);
  lower->SetRightMargin(0.04);

  upper->Draw();
  lower->Draw();

  upper->cd();

  mc.SetMarkerStyle(24);
  mc.SetMarkerSize(0.9);
  mc.SetMarkerColor(kBlack);
  mc.SetLineColor(kBlack);
  mc.SetLineWidth(2);

  if (dataGraph) {
    dataGraph->SetMarkerStyle(20);
    dataGraph->SetMarkerSize(1.0);
    dataGraph->SetMarkerColor(kRed + 1);
    dataGraph->SetLineColor(kRed + 1);
  }

  double maxVal = std::max(mc.GetMaximum(), dataHist.GetMaximum());
  if (maxVal <= 0.0)
    maxVal = 1.0;

  mc.SetMinimum(0.0);
  mc.SetMaximum(maxVal * 1.25);
  mc.GetXaxis()->SetLabelSize(0.0);
  mc.GetXaxis()->SetTitleSize(0.0);
  mc.GetYaxis()->SetTitleOffset(1.35);
  mc.Draw("E1");

  if (dataGraph)
    dataGraph->Draw("P SAME");
  else {
    dataHist.SetMarkerStyle(20);
    dataHist.SetMarkerColor(kRed + 1);
    dataHist.SetLineColor(kRed + 1);
    dataHist.Draw("E1 SAME");
  }

  auto* legend = new TLegend(0.53, 0.70, 0.92, 0.89);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(&mc, "Ditto", "lep");
  if (dataGraph)
    legend->AddEntry(dataGraph, "ALICE, pp #sqrt{s}=13 TeV", "lep");
  else
    legend->AddEntry(&dataHist, "ALICE, pp #sqrt{s}=13 TeV", "lep");
  legend->AddEntry((TObject*)nullptr, "primary charged particles, INEL", "");
  legend->Draw();

  lower->cd();
  ratio.SetMarkerStyle(24);
  ratio.SetMarkerSize(0.8);
  ratio.GetYaxis()->SetRangeUser(0.5, 1.5);
  ratio.GetYaxis()->SetNdivisions(505);
  ratio.GetYaxis()->SetTitleSize(0.10);
  ratio.GetYaxis()->SetTitleOffset(0.55);
  ratio.GetYaxis()->SetLabelSize(0.09);
  ratio.GetXaxis()->SetTitleSize(0.12);
  ratio.GetXaxis()->SetLabelSize(0.10);
  ratio.Draw("E1");

  auto* unity = new TLine(ratio.GetXaxis()->GetXmin(), 1.0, ratio.GetXaxis()->GetXmax(), 1.0);
  unity->SetLineStyle(2);
  unity->Draw("SAME");

  canvas->cd();
  canvas->Modified();
  canvas->Update();
  std::string printName = pdfName;
  if (first)
    printName += "(";
  else if (last)
    printName += ")";
  canvas->Print(printName.c_str());
}

} // namespace DittoValidation

void ValidateDitto(const char* inputFile = "Generation/Ditto.root",
                   const char* hepDataDir = "Validation/hepdata",
                   const char* outputFile = "Validation/Ditto_validation.root",
                   const char* outputPdf = "Validation/Ditto_vs_ALICE_13TeV.pdf",
                   bool downloadHEPData = true)
{
  using namespace DittoValidation;

  gStyle->SetOptStat(0);

  // --------------------------------------------------------------------------
  // Load HEPData references and create MC histograms with IDENTICAL binning.
  // --------------------------------------------------------------------------

  struct WorkingSet {
    const ReferenceSpec* ref = nullptr;
    HepDataObjects hepData;
    TH1D* mc = nullptr;
    TH1D* ratio = nullptr;
  };

  std::vector<WorkingSet> sets;
  sets.reserve(std::size(kReferences));

  for (const auto& ref : kReferences) {
    if (downloadHEPData)
      downloadTable(hepDataDir, ref);

    const std::string fileName = hepDataFileName(hepDataDir, ref);
    auto hd = openHEPData(fileName, ref);
    auto mc = makeMCHistogram(*hd.hist, ref);

    WorkingSet set;
    set.ref = &ref;
    set.hepData = std::move(hd);
    set.mc = std::move(mc);
    sets.push_back(std::move(set));
  }

  // --------------------------------------------------------------------------
  // Charged-particle dN_ch/deta reference: ALICE pp 13 TeV, INEL.
  // HEPData 70847, Table 1, y1.
  // --------------------------------------------------------------------------

  if (downloadHEPData)
    downloadChargedEtaTable(hepDataDir);

  auto chargedEtaHD =
    openChargedEtaHEPData(chargedEtaHepDataFileName(hepDataDir));
  auto* hMC_Charged_Eta = makeMCEtaHistogram(*chargedEtaHD.hist);

  // --------------------------------------------------------------------------
  // Input Ditto tree.
  // --------------------------------------------------------------------------

  std::unique_ptr<TFile> input(TFile::Open(inputFile, "READ"));
  if (!input || input->IsZombie())
    throw std::runtime_error(std::string("Cannot open Ditto file ") + inputFile);

  auto* tree = dynamic_cast<TTree*>(input->Get("T"));
  if (!tree)
    throw std::runtime_error("Could not find TTree 'T'");

  TClonesArray* particles = nullptr;
  if (tree->SetBranchAddress("Particles", &particles) < 0)
    throw std::runtime_error("Could not connect branch 'Particles'");

  const Long64_t nEvents = tree->GetEntries();
  std::cout << "Reading " << nEvents << " Ditto events\n";

  // --------------------------------------------------------------------------
  // One event loop fills all species and charged particle distributions.
  // --------------------------------------------------------------------------

  for (Long64_t iev = 0; iev < nEvents; ++iev) {
    tree->GetEntry(iev);

    if (!particles)
      continue;

    const int nParticles = particles->GetEntriesFast();
    for (int ip = 0; ip < nParticles; ++ip) {
      const auto* particle =
        static_cast<const TParticle*>(particles->UncheckedAt(ip));
      if (!particle)
        continue;

      // Current Ditto output uses positive status for final generator-level
      // particles.  Long-lived particles such as K0S/Lambda are not decayed here;
      // that is left to transport, exactly as in the generator workflow.
      if (particle->GetStatusCode() <= 0)
        continue;

      const int pdgCode = particle->GetPdgCode();

      // All charged particles eta distribution
      if (isCharged(pdgCode)) {
        const double eta = particle->Eta();
        if (std::isfinite(eta))
          hMC_Charged_Eta->Fill(eta);
      }

      const double y = rapidity(*particle);
      if (!std::isfinite(y) || std::abs(y) >= kYMax)
        continue;

      const int absPdg = std::abs(pdgCode);
      const double pt = particle->Pt();

      for (auto& set : sets) {
        if (absPdg == set.ref->absPdg)
          set.mc->Fill(pt);
      }
    }
  }

  // --------------------------------------------------------------------------
  // Normalize histograms.
  // --------------------------------------------------------------------------

  normalizeEtaDistribution(*hMC_Charged_Eta, nEvents);
  auto* hRatio_Charged_Eta =
    makeEtaRatio(*hMC_Charged_Eta, *chargedEtaHD.hist);

  for (auto& set : sets) {
    normalizeSpectrum(*set.mc, nEvents);
    set.ratio = makeRatio(*set.mc, *set.hepData.hist, *set.ref);
  }

  // --------------------------------------------------------------------------
  // Save ROOT objects and make comparison PDF.
  // --------------------------------------------------------------------------

  std::unique_ptr<TFile> output(TFile::Open(outputFile, "RECREATE"));
  if (!output || output->IsZombie())
    throw std::runtime_error(std::string("Cannot create output file ") + outputFile);

  for (std::size_t i = 0; i < sets.size(); ++i) {
    auto& set = sets[i];

    output->cd();
    auto* dir = output->mkdir(set.ref->shortName);
    dir->cd();

    set.mc->Write("Ditto");
    set.ratio->Write("Ratio_MC_Data");

    auto* dataClone = dynamic_cast<TH1*>(set.hepData.hist->Clone("HEPData_Central"));
    if (dataClone) {
      dataClone->SetDirectory(dir);
      dataClone->Write();
    }

    if (set.hepData.graph) {
      auto* graphClone = dynamic_cast<TGraphAsymmErrors*>(
        set.hepData.graph->Clone("HEPData_Graph"));
      if (graphClone)
        graphClone->Write();
    }

    drawComparison(*set.ref,
                   *set.mc,
                   *set.hepData.hist,
                   set.hepData.graph,
                   *set.ratio,
                   outputPdf,
                   i == 0,
                   false);
  }

  output->cd();
  auto* dirCharged = output->mkdir("Charged");
  dirCharged->cd();
  hMC_Charged_Eta->Write("Ditto_Eta");
  hRatio_Charged_Eta->Write("Ratio_MC_Data");

  auto* chargedDataClone =
    dynamic_cast<TH1*>(chargedEtaHD.hist->Clone("HEPData_Central"));
  if (chargedDataClone) {
    chargedDataClone->SetDirectory(dirCharged);
    chargedDataClone->Write();
  }

  if (chargedEtaHD.graph) {
    auto* chargedGraphClone = dynamic_cast<TGraphAsymmErrors*>(
      chargedEtaHD.graph->Clone("HEPData_Graph"));
    if (chargedGraphClone)
      chargedGraphClone->Write();
  }

  drawChargedEta(*hMC_Charged_Eta,
                 *chargedEtaHD.hist,
                 chargedEtaHD.graph,
                 *hRatio_Charged_Eta,
                 outputPdf,
                 sets.empty(),
                 true);

  output->Close();

  std::cout << "\nValidation complete.\n"
            << "  ROOT output : " << outputFile << "\n"
            << "  plots       : " << outputPdf << "\n"
            << "  identified : HEPData " << kHepDataRecord
            << " v" << kHepDataVersion << " (Tables 1, 2, 6)\n"
            << "  charged eta: HEPData " << kChargedEtaHepDataRecord
            << " v" << kChargedEtaHepDataVersion
            << " (Table 1, INEL / y1)\n"
            << "  selection   : identified hadrons |y| < " << kYMax << "\n"
            << "  observables : (1/N_INEL) d^2N/(dpT dy), (1/N_INEL) dN_ch/deta\n"
            << "\nNote: the HEPData tables exclude the common 2.6% normalization uncertainty.\n";
}
