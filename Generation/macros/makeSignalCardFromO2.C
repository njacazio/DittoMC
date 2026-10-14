/// Example signal-card writer for Ditto from O2 MC output.
///
/// Injects signal particles with their natural, Nch-dependent number per event.
/// All inputs are read from the AnalysisResults.root of alice3-qa-singleparticle.
/// Every top-level directory containing AnalysedParticlePdgCode is one task
/// instance (e.g. one per --workflow-suffix) and becomes one signal species:
///   <directory>/AnalysedParticlePdgCode  PDG code the instance was run with
///   <directory>/AnalysedParticleMass     mass of that particle
///   <directory>/particle/NSignalVsNch    P(N_signal | Nch), summed over the Nch
///                                        bins of each tune activity class
///   <directory>/particle/Pt              pT shape
///   <directory>/particle/Eta             eta shape
/// The pT and eta shapes are used identically in every activity class.
///
/// With addAntiParticles, the antiparticle (-PDG) of every species is injected
/// as an additional, independent species with the same multiplicity, pT and eta
/// distributions. The task only counts the exact pdgCode, so this assumes
/// particle/antiparticle symmetry. Species whose antiparticle is unknown to
/// TDatabasePDG (self-conjugate, e.g. the J/psi) are not duplicated.

#include "DittoSignal.h"
#include "DittoTune.h"

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TParameter.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

template <typename T>
T* loadHistogram(TFile& file, const std::string& path, const std::string& newName)
{
  auto* hist = file.Get<T>(path.c_str());
  if (!hist) {
    throw std::runtime_error("Could not read " + path + " from " + file.GetName());
  }
  hist->SetDirectory(nullptr);
  hist->SetName(newName.c_str());
  return hist;
}

std::vector<double> getBinEdges(const TAxis& axis)
{
  std::vector<double> edges;
  for (int iBin = 1; iBin <= axis.GetNbins() + 1; ++iBin) {
    edges.push_back(axis.GetBinLowEdge(iBin));
  }
  return edges;
}

// Fills P(N_signal | activity class) for one class by summing the N_signal
// distributions of the input Nch bins whose centres lie in [nchMin, nchMax).
// The card must have the same N_signal binning as the input. A class without
// any input events is not injected (N_signal = 0).
void fillSignalMultiplicity(TH2D& card, int iActivity, const TH2& nSignalVsNch, double nchMin, double nchMax)
{
  double nEvents = 0.0;
  double nOverflowEvents = 0.0;
  for (int iNch = 1; iNch <= nSignalVsNch.GetNbinsX(); ++iNch) {
    const double nch = nSignalVsNch.GetXaxis()->GetBinCenter(iNch);
    if (nch < nchMin || nch >= nchMax) {
      continue;
    }
    for (int iSignal = 1; iSignal <= nSignalVsNch.GetNbinsY(); ++iSignal) {
      const double content = nSignalVsNch.GetBinContent(iNch, iSignal);
      card.SetBinContent(iActivity, iSignal, card.GetBinContent(iActivity, iSignal) + content);
      nEvents += content;
    }
    nOverflowEvents += nSignalVsNch.GetBinContent(iNch, nSignalVsNch.GetNbinsY() + 1);
  }

  if (nOverflowEvents > 0.0) {
    Printf("Warning: %g events with N_signal above the %s axis dropped for Nch [%g, %g)", nOverflowEvents, nSignalVsNch.GetName(), nchMin, nchMax);
  }
  if (nEvents <= 0.0) {
    Printf("Warning: no events in %s for Nch [%g, %g), no injection in this class", nSignalVsNch.GetName(), nchMin, nchMax);
    card.SetBinContent(iActivity, card.GetYaxis()->FindBin(0.0), 1.0);
  }
}

// The task stores single values (PDG code, mass) as the content of a one-bin
// histogram, set once per job. Merged outputs (hadd, Hyperloop) sum the bin over
// jobs, so divide by the number of entries, which SetBinContent increments once
// per job.
double getStoredValue(TFile& file, const std::string& path)
{
  std::unique_ptr<TH1> hist(loadHistogram<TH1>(file, path, "storedValue"));
  return hist->GetBinContent(1) / std::max(1.0, hist->GetEntries());
}

std::pair<int, double> getAnalysedParticle(TFile& file, const std::string& directory)
{
  const double pdg = getStoredValue(file, directory + "/AnalysedParticlePdgCode");
  const double mass = getStoredValue(file, directory + "/AnalysedParticleMass");
  if (pdg == 0.0 || std::abs(pdg - std::round(pdg)) > 1e-6) {
    throw std::runtime_error("Invalid PDG code in " + directory + "/AnalysedParticlePdgCode: " + std::to_string(pdg));
  }
  if (!(mass > 0.0)) {
    throw std::runtime_error("Invalid mass in " + directory + "/AnalysedParticleMass: " + std::to_string(mass) + " (PDG code unknown to the task?)");
  }
  return {static_cast<int>(std::lround(pdg)), mass};
}

// One analysed species, read from an alice3-qa-singleparticle task directory.
struct SignalSpecies {
  std::string directory; // alice3-qa-singleparticle instance
  int pdg;               // from <directory>/AnalysedParticlePdgCode
  double mass;           // GeV/c^2, from <directory>/AnalysedParticleMass
  TH1* ptShape;
  TH1* etaShape;
  TH2* nSignalVsNch; // events vs (Nch, N_signal)

  SignalSpecies(TFile& file, const std::string& directoryName)
    : directory(directoryName)
  {
    std::tie(pdg, mass) = getAnalysedParticle(file, directory);
    const std::string tag = std::to_string(pdg);
    Printf("%s: PDG %d, mass %g GeV/c^2", directory.c_str(), pdg, mass);
    ptShape = loadHistogram<TH1>(file, directory + "/particle/Pt", "pt_" + tag);
    etaShape = loadHistogram<TH1>(file, directory + "/particle/Eta", "eta_" + tag);
    nSignalVsNch = loadHistogram<TH2>(file, directory + "/particle/NSignalVsNch", "nSignalVsNch_" + tag);
  }
};

// One species per alice3-qa-singleparticle directory in the file, i.e. per
// top-level directory that contains AnalysedParticlePdgCode.
std::vector<SignalSpecies> findSignalSpecies(TFile& file)
{
  std::vector<SignalSpecies> species;
  for (auto* key : TRangeDynCast<TKey>(file.GetListOfKeys())) {
    auto* directory = file.GetDirectory(key->GetName());
    if (!directory || !directory->FindKey("AnalysedParticlePdgCode")) {
      continue;
    }
    species.emplace_back(file, key->GetName());
  }
  if (species.empty()) {
    throw std::runtime_error(std::string("No alice3-qa-singleparticle directories found in ") + file.GetName());
  }
  return species;
}

void makeSignalCardFromO2(const char* tuneFile = "https://github.com/njacazio/DittoMC/releases/download/v1.0.1/Ditto_tune_pythia8_inel_136tev.root",
                          const char* analysisResultsFile = "AnalysisResults.root",
                          const char* outputFile = "signal_injected.root")
{
  const int enableDecay = 0;
  const bool addAntiParticles = true;

  std::unique_ptr<TFile> analysisResults(TFile::Open(analysisResultsFile));
  if (!analysisResults || analysisResults->IsZombie()) {
    throw std::runtime_error(std::string("Could not open ") + analysisResultsFile);
  }
  const auto species = findSignalSpecies(*analysisResults);
  analysisResults->Close();

  // Antiparticles share the histograms of their particle
  std::vector<SignalSpecies> injectedSpecies;
  for (const auto& s : species) {
    injectedSpecies.push_back(s);
    if (!addAntiParticles) {
      continue;
    }
    if (!TDatabasePDG::Instance()->GetParticle(-s.pdg)) {
      Printf("%s: PDG %d has no antiparticle in TDatabasePDG (self-conjugate?), not adding one", s.directory.c_str(), s.pdg);
      continue;
    }
    Printf("%s: adding antiparticle PDG %d", s.directory.c_str(), -s.pdg);
    auto antiparticle = s;
    antiparticle.pdg = -s.pdg;
    injectedSpecies.push_back(antiparticle);
  }

  // Ditto requires every PDG to appear once
  for (std::size_t i = 0; i < injectedSpecies.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      if (injectedSpecies[i].pdg == injectedSpecies[j].pdg) {
        throw std::runtime_error("PDG " + std::to_string(injectedSpecies[i].pdg) + " from both " + injectedSpecies[j].directory + " and " + injectedSpecies[i].directory);
      }
    }
  }

  auto tune = Ditto::Tune::load(tuneFile);
  if (!tune || tune->mActivityEdges.size() < 2) {
    throw std::runtime_error("Could not load valid Ditto activity edges");
  }
  const auto& activityEdges = tune->mActivityEdges;
  const int nActivityBins = static_cast<int>(activityEdges.size()) - 1;

  std::vector<int> pdgs;
  std::vector<double> masses;
  for (const auto& s : injectedSpecies) {
    pdgs.push_back(s.pdg);
    masses.push_back(s.mass);
  }

  TFile* output = TFile::Open(outputFile, "RECREATE");
  if (!output || output->IsZombie()) {
    throw std::runtime_error(std::string("Could not create ") + outputFile);
  }

  TH1D* massHist = Ditto::SignalInjector::makeSignalMassHistogram(*tune, pdgs, masses);
  massHist->Write();

  TParameter<int> signalDecay("signalDecay", enableDecay);
  signalDecay.Write();

  for (const auto& s : injectedSpecies) {
    // P(N_signal | activity, PDG), same N_signal binning as the input
    TH2D* multiplicityHist = Ditto::SignalInjector::makeSignalMultiplicityHistogram(*tune, s.pdg, getBinEdges(*s.nSignalVsNch->GetYaxis()));
    // P(pT | activity, PDG)
    TH2D* ptHist = Ditto::SignalInjector::makePtHistogram(*tune, s.pdg, getBinEdges(*s.ptShape->GetXaxis()));
    // P(eta | activity, PDG)
    TH2D* etaHist = Ditto::SignalInjector::makeEtaHistogram(*tune, s.pdg, getBinEdges(*s.etaShape->GetXaxis()));

    for (int iActivity = 1; iActivity <= nActivityBins; ++iActivity) {
      fillSignalMultiplicity(*multiplicityHist, iActivity, *s.nSignalVsNch, activityEdges[iActivity - 1], activityEdges[iActivity]);

      for (int iPt = 1; iPt <= ptHist->GetNbinsY(); ++iPt) {
        ptHist->SetBinContent(iActivity, iPt, s.ptShape->GetBinContent(iPt));
      }
      for (int iEta = 1; iEta <= etaHist->GetNbinsY(); ++iEta) {
        etaHist->SetBinContent(iActivity, iEta, s.etaShape->GetBinContent(iEta));
      }
    }

    multiplicityHist->Write();
    ptHist->Write();
    etaHist->Write();
  }

  output->Close();

  Printf("Wrote signal card to %s", outputFile);
}
