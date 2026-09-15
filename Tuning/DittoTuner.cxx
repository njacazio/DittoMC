///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTuner.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Implementation of the Ditto teacher tuner.
///

#include "DittoTuner.h"

#include "DittoTune.h"

// ACLiC/rootcling only needs to see the public declarations from the header.
// Hide all implementation details from dictionary generation.
//
// This keeps PYTHIA and all runtime-only implementation details out of the
// ROOT dictionary while the persistent Tune object remains fully streamable.
#ifndef __ROOTCLING__

#include <Pythia8/Pythia.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Ditto
{

namespace
{

std::string cardStem(const std::string& card)
{
  const std::size_t slash = card.find_last_of("/\\");
  const std::string filename = slash == std::string::npos ? card : card.substr(slash + 1);
  const std::size_t dot = filename.find_last_of('.');
  return dot == std::string::npos ? filename : filename.substr(0, dot);
}

std::string cardPath(const std::string& card)
{
  const std::size_t slash = card.find_last_of("/\\");
  return slash == std::string::npos ? "" : card.substr(0, slash + 1);
}

std::string readTextFile(const std::string& fileName)
{
  std::ifstream input(fileName);

  if (!input) {
    throw std::runtime_error("Ditto::Tuner: could not read PYTHIA card: " + fileName);
  }

  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

std::uint64_t compositionKey(int nch, int nSelected)
{
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(nch)) << 32) | static_cast<std::uint32_t>(nSelected);
}

int compositionKeyNch(std::uint64_t key)
{
  return static_cast<int>(static_cast<std::uint32_t>(key >> 32));
}

int compositionKeyNSelected(std::uint64_t key)
{
  return static_cast<int>(static_cast<std::uint32_t>(key & 0xffffffffULL));
}

struct CompositionTemplate {
  std::vector<unsigned short> centralCounts;
  std::vector<unsigned short> otherCounts;
};

struct CompositionReservoir {
  std::uint64_t eventsSeen = 0;
  std::vector<CompositionTemplate> templates;
};

} // namespace

TunerConfig::TunerConfig()
{
  species = {22,
             11,
             -11,
             13,
             -13,
             211,
             -211,
             111,
             321,
             -321,
             310,
             130,
             2212,
             -2212,
             2112,
             -2112,
             3122,
             -3122};
}

struct TunerImpl {
  explicit TunerImpl(const TunerConfig& cfg)
    : config(cfg),
      compositionRng(cfg.compositionReservoirSeed)
  {
  }

  TunerConfig config;

  std::string outputFile;

  std::unique_ptr<Pythia8::Pythia> pythia;

  Tune tune;
  std::unordered_map<int, std::size_t> speciesIndex;

  std::mt19937_64 compositionRng;
  std::unordered_map<std::uint64_t, CompositionReservoir> compositionReservoirs;

  std::uint64_t generatedEvents = 0;
  std::uint64_t attemptedEvents = 0;
  std::uint64_t activityOverflowEvents = 0;
  std::uint64_t ptOverflowParticles = 0;
  std::uint64_t speciesMultiplicityOverflowEvents = 0;

  bool ran = false;

  void validateConfig() const
  {
    if (config.pythiaCard.empty()) {
      throw std::invalid_argument("Ditto::Tuner: pythiaCard must not be empty");
    }

    if (config.nEvents == 0) {
      throw std::invalid_argument("Ditto::Tuner: nEvents must be > 0");
    }

    if (config.activityEtaMax <= 0.0) {
      throw std::invalid_argument("Ditto::Tuner: activityEtaMax must be > 0");
    }

    if (config.particleEtaMax <= 0.0) {
      throw std::invalid_argument("Ditto::Tuner: particleEtaMax must be > 0");
    }

    if (config.finalStatus <= 0) {
      throw std::invalid_argument("Ditto::Tuner: finalStatus must be positive");
    }

    if (config.activityEdges.size() < 2) {
      throw std::invalid_argument("Ditto::Tuner: at least two activity edges are required");
    }

    if (!std::is_sorted(config.activityEdges.begin(), config.activityEdges.end())) {
      throw std::invalid_argument("Ditto::Tuner: activityEdges must be sorted");
    }

    for (std::size_t i = 1; i < config.activityEdges.size(); ++i) {
      if (config.activityEdges[i] <= config.activityEdges[i - 1]) {
        throw std::invalid_argument("Ditto::Tuner: activityEdges must be strictly increasing");
      }
    }

    if (config.selectedMultiplicityEdges.size() < 2) {
      throw std::invalid_argument("Ditto::Tuner: at least two selectedMultiplicityEdges are required");
    }

    for (std::size_t i = 1; i < config.selectedMultiplicityEdges.size(); ++i) {
      if (config.selectedMultiplicityEdges[i] <= config.selectedMultiplicityEdges[i - 1]) {
        throw std::invalid_argument("Ditto::Tuner: selectedMultiplicityEdges must be strictly increasing");
      }
    }

    if (config.species.empty()) {
      throw std::invalid_argument("Ditto::Tuner: species list must not be empty");
    }

    if (config.maxNch <= 0 || config.maxSelectedMultiplicity <= 0 || config.maxSpeciesMultiplicity <= 0) {
      throw std::invalid_argument("Ditto::Tuner: multiplicity histogram maxima must be > 0");
    }

    if (config.nPtBins <= 0 || config.ptMax <= 0.0) {
      throw std::invalid_argument("Ditto::Tuner: invalid pT binning");
    }

    if (config.nEtaBins <= 0) {
      throw std::invalid_argument("Ditto::Tuner: nEtaBins must be > 0");
    }

    if (config.maxCompositionTemplatesPerPair == 0) {
      throw std::invalid_argument("Ditto::Tuner: maxCompositionTemplatesPerPair must be > 0");
    }

    if (config.maxAttemptsFactor < 1.0) {
      throw std::invalid_argument("Ditto::Tuner: maxAttemptsFactor must be >= 1");
    }
  }

  void configurePythia()
  {
    if (!pythia->readFile(config.pythiaCard)) {
      throw std::runtime_error("Ditto::Tuner: could not read PYTHIA card: " + config.pythiaCard);
    }

    const int frameType = pythia->settings.mode("Beams:frameType");

    if (frameType != 1) {
      throw std::invalid_argument(
        "Ditto::Tuner: only Beams:frameType = 1 is currently supported");
    }

    const double sqrtSNN = pythia->settings.parm("Beams:eCM");

    if (!std::isfinite(sqrtSNN) ||
        sqrtSNN <= 0.0) {
      throw std::invalid_argument("Ditto::Tuner: Beams:eCM must be > 0");
    }

    if (!pythia->init()) {
      throw std::runtime_error("Ditto::Tuner: PYTHIA initialization failed");
    }

    for (const int pdg : config.species) {
      if (!pythia->particleData.isParticle(pdg)) {
        throw std::invalid_argument("Ditto::Tuner: unknown PDG code " + std::to_string(pdg));
      }
    }
  }

  void bookTune()
  {
    tune.initialize(config.activityEtaMax,
                    config.particleEtaMax,
                    config.activityEdges,
                    config.selectedMultiplicityEdges,
                    config.maxNch,
                    config.maxSelectedMultiplicity,
                    config.maxSpeciesMultiplicity,
                    config.nPtBins,
                    config.ptMax,
                    config.nEtaBins,
                    config.species);

    tune.teacher = "PYTHIA8";
    tune.pythiaCard = config.pythiaCard;
    tune.pythiaCardContent = readTextFile(config.pythiaCard);

    tune.beamIdA = pythia->settings.mode("Beams:idA");

    tune.beamIdB = pythia->settings.mode("Beams:idB");

    tune.beamFrameType = pythia->settings.mode("Beams:frameType");

    tune.sqrtSNN = pythia->settings.parm("Beams:eCM");

    tune.azimuthModel = "uniform";
    tune.finalStatus = config.finalStatus;

    speciesIndex.clear();
    speciesIndex.reserve(config.species.size());

    for (std::size_t i = 0; i < config.species.size(); ++i) {
      const int pdg = config.species[i];

      speciesIndex[pdg] = i;

      auto* entry = tune.speciesAt(static_cast<int>(i));

      entry->particleName = pythia->particleData.name(pdg);

      entry->mass = pythia->particleData.m0(pdg);

      entry->chargeType = pythia->particleData.chargeType(pdg);
    }
  }

  int activityClass(double nch) const
  {
    return tune.activityClass(nch);
  }

  int selectedMultiplicityClass(double nSelected) const
  {
    return tune.selectedMultiplicityClass(nSelected);
  }

  void storeCompositionTemplate(int nch,
                                int nSelected,
                                const std::vector<int>& centralCounts,
                                const std::vector<int>& otherCounts)
  {
    const std::uint64_t key = compositionKey(nch, nSelected);

    auto& reservoir = compositionReservoirs[key];

    ++reservoir.eventsSeen;

    std::size_t destination = 0;

    if (reservoir.templates.size() < config.maxCompositionTemplatesPerPair) {
      destination = reservoir.templates.size();
      reservoir.templates.emplace_back();
    } else {
      std::uniform_int_distribution<std::uint64_t> distribution(0, reservoir.eventsSeen - 1);

      const std::uint64_t selected = distribution(compositionRng);

      if (selected >= config.maxCompositionTemplatesPerPair) {
        return;
      }

      destination = static_cast<std::size_t>(selected);
    }

    auto& composition = reservoir.templates[destination];

    composition.centralCounts.resize(centralCounts.size());
    composition.otherCounts.resize(otherCounts.size());

    for (std::size_t i = 0; i < centralCounts.size(); ++i) {
      if (centralCounts[i] < 0 ||
          otherCounts[i] < 0 ||
          centralCounts[i] > static_cast<int>(std::numeric_limits<unsigned short>::max()) ||
          otherCounts[i] > static_cast<int>(std::numeric_limits<unsigned short>::max())) {
        throw std::runtime_error("Ditto::Tuner: composition multiplicity exceeds unsigned-short storage");
      }

      composition.centralCounts[i] = static_cast<unsigned short>(centralCounts[i]);

      composition.otherCounts[i] = static_cast<unsigned short>(otherCounts[i]);
    }
  }

  void serializeCompositionTemplates()
  {
    tune.compositionTemplateCapPerPair = static_cast<std::uint64_t>(config.maxCompositionTemplatesPerPair);

    tune.compositionPairNch.clear();
    tune.compositionPairNSelected.clear();
    tune.compositionPairOffsets.clear();
    tune.compositionPairEventsSeen.clear();
    tune.compositionCentralCounts.clear();
    tune.compositionOtherCounts.clear();

    std::vector<std::uint64_t> keys;
    keys.reserve(compositionReservoirs.size());

    for (const auto& item : compositionReservoirs) {
      if (!item.second.templates.empty()) {
        keys.push_back(item.first);
      }
    }

    std::sort(keys.begin(), keys.end());

    tune.compositionPairNch.reserve(keys.size());
    tune.compositionPairNSelected.reserve(keys.size());
    tune.compositionPairEventsSeen.reserve(keys.size());
    tune.compositionPairOffsets.reserve(keys.size() + 1);

    tune.compositionPairOffsets.push_back(0);

    const std::size_t nSpecies = config.species.size();

    std::uint64_t templateOffset = 0;

    for (const std::uint64_t key : keys) {
      const auto found = compositionReservoirs.find(key);

      if (found == compositionReservoirs.end()) {
        continue;
      }

      const auto& reservoir = found->second;

      tune.compositionPairNch.push_back(compositionKeyNch(key));

      tune.compositionPairNSelected.push_back(compositionKeyNSelected(key));

      tune.compositionPairEventsSeen.push_back(reservoir.eventsSeen);

      for (const auto& composition : reservoir.templates) {
        if (composition.centralCounts.size() != nSpecies || composition.otherCounts.size() != nSpecies) {
          throw std::runtime_error("Ditto::Tuner: inconsistent composition-template species count");
        }

        tune.compositionCentralCounts.insert(tune.compositionCentralCounts.end(),
                                             composition.centralCounts.begin(),
                                             composition.centralCounts.end());

        tune.compositionOtherCounts.insert(tune.compositionOtherCounts.end(),
                                           composition.otherCounts.begin(),
                                           composition.otherCounts.end());
      }

      templateOffset += static_cast<std::uint64_t>(reservoir.templates.size());

      tune.compositionPairOffsets.push_back(templateOffset);
    }
  }

  void analyzeCurrentEvent()
  {
    const auto& event = pythia->event;

    int nch = 0;

    // ----------------------------------------------------------------------
    // First pass: define the event activity.
    // ----------------------------------------------------------------------

    for (int i = 1; i < event.size(); ++i) {
      const auto& p = event[i];

      if (!p.isFinal())
        continue;

      const double eta = p.eta();

      if (!std::isfinite(eta) || std::abs(eta) >= config.activityEtaMax) {
        continue;
      }

      if (p.isCharged())
        ++nch;
    }

    tune.hNch.Fill(nch);

    if (activityClass(nch) < 0) {
      ++activityOverflowEvents;
      return;
    }

    std::vector<int> speciesCounts(config.species.size(), 0);
    std::vector<int> centralChargedSpeciesCounts(config.species.size(), 0);
    std::vector<int> otherSpeciesCounts(config.species.size(), 0);

    int nSelected = 0;
    int nSelectedCentralCharged = 0;

    // ----------------------------------------------------------------------
    // Second pass: composition and single-particle kinematics conditioned
    // on Nch.
    // ----------------------------------------------------------------------

    for (int i = 1; i < event.size(); ++i) {
      const auto& p = event[i];

      if (!p.isFinal())
        continue;

      const auto found = speciesIndex.find(p.id());

      if (found == speciesIndex.end())
        continue;

      const double eta = p.eta();

      if (!std::isfinite(eta) || std::abs(eta) >= config.particleEtaMax) {
        continue;
      }

      const std::size_t iSpecies = found->second;

      ++speciesCounts[iSpecies];
      ++nSelected;

      const bool centralCharged = p.isCharged() && std::abs(eta) < config.activityEtaMax;

      if (centralCharged) {
        ++centralChargedSpeciesCounts[iSpecies];
        ++nSelectedCentralCharged;
      } else {
        ++otherSpeciesCounts[iSpecies];
      }

      const double pt = p.pT();

      if (pt >= config.ptMax)
        ++ptOverflowParticles;

      auto* entry = tune.speciesAt(static_cast<int>(iSpecies));

      entry->hPtVsActivity.Fill(nch, pt);
      entry->hEtaVsActivity.Fill(nch, eta);

      if (centralCharged) {
        entry->hCentralChargedPtSumVsNch.Fill(nch, pt);
        entry->hCentralChargedPtCountVsNch.Fill(nch);
      } else {
        entry->hOtherPtSumVsNch.Fill(nch, pt);
        entry->hOtherPtCountVsNch.Fill(nch);
      }
    }

    tune.hNSelected.Fill(nSelected);
    tune.hNchSelected.Fill(nSelectedCentralCharged);

    const double selectedMultiplicityPoint[2] = {static_cast<double>(nch),
                                                 static_cast<double>(nSelected)};
    const double selectedChargedPoint[2] = {static_cast<double>(nch),
                                            static_cast<double>(nSelectedCentralCharged)};

    tune.hNSelectedVsNch.Fill(selectedMultiplicityPoint);
    tune.hNchSelectedVsNch.Fill(selectedChargedPoint);

    const int selectedClass = selectedMultiplicityClass(nSelected);

    if (selectedClass < 0) {
      ++tune.selectedMultiplicityOverflowEvents;
    } else {
      tune.hEventsVsActivitySelected.Fill(nch, nSelected);
    }

    if (nSelectedCentralCharged != nch) {
      ++tune.centralChargedCoverageMismatchEvents;

      if (nSelectedCentralCharged < nch) {
        tune.centralChargedCoverageMissingParticles += static_cast<std::uint64_t>(nch - nSelectedCentralCharged);
      }
    }

    // Fill once per event for every species, including N_species = 0.
    //
    // hCountVsActivity retains the original per-species multiplicity model
    // for diagnostics and backwards comparison.
    //
    // The TH2D composition tables retain the mean species fractions used by
    // tune format v2. The sparse two-dimensional tables retain the complete
    // central / other multiplicity PMFs versus the Nch activity class. They are
    // filled once per event, including zero, so that they represent
    //
    //   P(N_species | activity class).
    for (std::size_t i = 0; i < speciesCounts.size(); ++i) {
      const int count = speciesCounts[i];

      if (count > config.maxSpeciesMultiplicity)
        ++speciesMultiplicityOverflowEvents;

      auto* entry = tune.speciesAt(static_cast<int>(i));

      entry->hCountVsActivity.Fill(nch, count);

      const int centralCount = centralChargedSpeciesCounts[i];
      const int otherCount = otherSpeciesCounts[i];

      const double centralMultiplicity[2] = {static_cast<double>(nch),
                                             static_cast<double>(centralCount)};
      const double otherMultiplicity[2] = {static_cast<double>(nch),
                                           static_cast<double>(otherCount)};

      entry->hCentralChargedMultiplicityVsActivity.Fill(centralMultiplicity);
      entry->hOtherMultiplicityVsActivity.Fill(otherMultiplicity);

      if (selectedClass >= 0) {
        if (centralCount > 0) {
          entry->hCentralChargedCountVsActivitySelected.Fill(nch, nSelected, centralCount);
        }

        if (otherCount > 0) {
          entry->hOtherCountVsActivitySelected.Fill(nch, nSelected, otherCount);
        }
      }
    }

    if (selectedClass >= 0) {
      storeCompositionTemplate(nch, nSelected, centralChargedSpeciesCounts, otherSpeciesCounts);
    }
  }

  void saveTune()
  {
    tune.nEvents = generatedEvents;
    tune.nAttempts = attemptedEvents;

    tune.activityOverflowEvents = activityOverflowEvents;
    tune.ptOverflowParticles = ptOverflowParticles;
    tune.speciesMultiplicityOverflowEvents = speciesMultiplicityOverflowEvents;

    serializeCompositionTemplates();

    tune.finalize();
    tune.save(outputFile);
  }

  void run()
  {
    if (ran) {
      throw std::runtime_error("Ditto::Tuner: run() may only be called once per Tuner instance");
    }

    ran = true;

    const std::uint64_t maxAttempts = static_cast<std::uint64_t>(std::ceil(config.maxAttemptsFactor * static_cast<double>(config.nEvents)));

    const auto start = std::chrono::steady_clock::now();

    while (generatedEvents < config.nEvents) {
      if (attemptedEvents >= maxAttempts) {
        throw std::runtime_error("Ditto::Tuner: too many failed PYTHIA attempts");
      }

      ++attemptedEvents;

      if (!pythia->next())
        continue;

      analyzeCurrentEvent();

      ++generatedEvents;

      if (config.progressEvery > 0 && generatedEvents % config.progressEvery == 0) {
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        const double eventsPerSecond = elapsed > 0.0 ? static_cast<double>(generatedEvents) / elapsed : 0.0;
        const double etaSeconds = eventsPerSecond > 0.0 ? static_cast<double>(config.nEvents - generatedEvents) / eventsPerSecond : 0.0;

        std::cout << "Ditto tuner: "
                  << generatedEvents
                  << " / "
                  << config.nEvents
                  << " successful PYTHIA events (ETA: " << etaSeconds << " s)\n";
      }
    }

    const auto end = std::chrono::steady_clock::now();

    const double elapsed = std::chrono::duration<double>(end - start).count();
    tune.generationTimeSeconds = elapsed;

    std::cout << "\nDitto tuner generation complete\n"
              << "  PYTHIA card      : "
              << config.pythiaCard
              << "\n"
              << "  beams            : "
              << tune.beamIdA
              << " + "
              << tune.beamIdB
              << "\n"
              << "  sqrt(s_NN)       : "
              << tune.sqrtSNN
              << " GeV\n"
              << "  events           : "
              << generatedEvents
              << "\n"
              << "  attempts         : "
              << attemptedEvents
              << "\n"
              << "  wall time        : "
              << std::fixed
              << std::setprecision(3)
              << elapsed
              << " s\n";

    if (elapsed > 0.0) {
      std::cout << "  teacher rate     : "
                << std::setprecision(2)
                << static_cast<double>(generatedEvents) / elapsed
                << " events/s\n";
    }

    if (activityOverflowEvents > 0) {
      std::cout << "  WARNING activity overflow events: "
                << activityOverflowEvents
                << " (extend activityEdges)\n";
    }

    if (tune.selectedMultiplicityOverflowEvents > 0) {
      std::cout << "  WARNING selected-multiplicity overflow events: "
                << tune.selectedMultiplicityOverflowEvents
                << " (extend selectedMultiplicityEdges)\n";
    }

    if (ptOverflowParticles > 0) {
      std::cout << "  WARNING pT overflow particles: "
                << ptOverflowParticles
                << " (increase ptMax if relevant)\n";
    }

    if (tune.centralChargedCoverageMismatchEvents > 0) {
      const double eventFraction = static_cast<double>(tune.centralChargedCoverageMismatchEvents) /
                                   static_cast<double>(generatedEvents);

      std::cout << "  WARNING selected species do not cover all central charged particles in "
                << tune.centralChargedCoverageMismatchEvents
                << " events ("
                << std::setprecision(4)
                << 100.0 * eventFraction
                << "%), missing "
                << tune.centralChargedCoverageMissingParticles
                << " charged particles in total\n";
    }

    if (speciesMultiplicityOverflowEvents > 0) {
      std::cout << "  WARNING species-count overflow fills: "
                << speciesMultiplicityOverflowEvents
                << " (increase maxSpeciesMultiplicity)\n";
    }

    if (config.printPythiaStatistics)
      pythia->stat();

    std::cout << "\nComputing conditional probabilities " << "and writing tune...\n";

    saveTune();

    std::cout << "  composition pairs     : "
              << tune.numberOfCompositionPairs()
              << "\n"
              << "  templates stored      : "
              << tune.numberOfCompositionTemplates()
              << "\n"
              << "  template cap / pair   : "
              << tune.compositionTemplateCapPerPair
              << "\n"
              << "Tune written to: "
              << outputFile
              << "\n";
  }
};

Tuner::Tuner(const TunerConfig& config)
{
  // Keep construction exception-safe even though the public class stores a
  // raw opaque pointer for ROOT/rootcling compatibility.
  auto impl = std::make_unique<TunerImpl>(config);

  impl->validateConfig();

  const std::string cardDir = cardPath(impl->config.pythiaCard);
  const std::string cardFileStem = cardStem(impl->config.pythiaCard);
  const std::string outputDir = (cardDir.empty() ? "." : cardDir) + "/../tunes/";
  if (!std::filesystem::exists(outputDir)) {
    std::filesystem::create_directories(outputDir);
  }
  impl->outputFile = outputDir + "Ditto_tune_" + cardFileStem + ".root";

  // Use PYTHIA's standard construction / installation lookup. In
  // particular, PYTHIA8DATA is honored by PYTHIA itself when it is set.
  impl->pythia = std::make_unique<Pythia8::Pythia>();

  impl->configurePythia();
  impl->bookTune();

  // Transfer ownership only after all potentially throwing initialization
  // has completed.
  mImpl = impl.release();
}

Tuner::~Tuner()
{
  delete mImpl;
  mImpl = nullptr;
}

void Tuner::run()
{
  if (!mImpl)
    throw std::runtime_error("Ditto::Tuner: invalid implementation");

  mImpl->run();
}

const TunerConfig& Tuner::config() const
{
  if (!mImpl)
    throw std::runtime_error("Ditto::Tuner: invalid implementation");

  return mImpl->config;
}

std::uint64_t Tuner::generatedEvents() const
{
  return mImpl ? mImpl->generatedEvents : 0;
}

std::uint64_t Tuner::attemptedEvents() const
{
  return mImpl ? mImpl->attemptedEvents : 0;
}

const Tune& Tuner::tune() const
{
  if (!mImpl)
    throw std::runtime_error("Ditto::Tuner: invalid implementation");

  return mImpl->tune;
}

} // namespace Ditto

#endif // __ROOTCLING__
