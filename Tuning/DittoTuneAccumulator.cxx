///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTuneAccumulator.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/19
/// \brief  Generator-independent accumulation of Ditto tuning observables.
///

#include "DittoTuneAccumulator.h"

#include "DittoTune.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Ditto
{

namespace
{

std::uint64_t compositionKey(int nch, int nSelected)
{
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(nch)) << 32) |
         static_cast<std::uint32_t>(nSelected);
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

struct TuneAccumulatorImpl {
  TuneAccumulatorImpl(const TuneAccumulatorConfig& cfg,
                      const std::vector<ParticleDefinition>& particleDefinitions)
    : config(cfg),
      compositionRng(cfg.mCompositionReservoirSeed)
  {
    validateConfig();
    bookTune(particleDefinitions);
  }

  TuneAccumulatorConfig config;
  Tune tune;
  std::unordered_map<int, std::size_t> speciesIndex;

  std::mt19937_64 compositionRng;
  std::unordered_map<std::uint64_t, CompositionReservoir> compositionReservoirs;

  std::uint64_t processedEvents = 0;
  std::uint64_t activityOverflowEvents = 0;
  std::uint64_t ptOverflowParticles = 0;
  std::uint64_t speciesMultiplicityOverflowEvents = 0;
  bool isFinalized = false;

  void validateConfig() const
  {
    if (config.mActivityEtaMax <= 0.0) {
      throw std::invalid_argument("Ditto::TuneAccumulator: mActivityEtaMax must be > 0");
    }
    if (config.mParticleEtaMax <= 0.0) {
      throw std::invalid_argument("Ditto::TuneAccumulator: particleEtaMax must be > 0");
    }
    if (config.mFinalStatus <= 0) {
      throw std::invalid_argument("Ditto::TuneAccumulator: finalStatus must be positive");
    }

    if (config.mActivityEdges.size() < 2) {
      throw std::invalid_argument("Ditto::TuneAccumulator: at least two activity edges are required");
    }
    for (std::size_t i = 1; i < config.mActivityEdges.size(); ++i) {
      if (config.mActivityEdges[i] <= config.mActivityEdges[i - 1]) {
        throw std::invalid_argument("Ditto::TuneAccumulator: activityEdges must be strictly increasing");
      }
    }

    if (config.mSelectedMultiplicityEdges.size() < 2) {
      throw std::invalid_argument(
        "Ditto::TuneAccumulator: at least two selectedMultiplicityEdges are required");
    }
    for (std::size_t i = 1; i < config.mSelectedMultiplicityEdges.size(); ++i) {
      if (config.mSelectedMultiplicityEdges[i] <= config.mSelectedMultiplicityEdges[i - 1]) {
        throw std::invalid_argument(
          "Ditto::TuneAccumulator: selectedMultiplicityEdges must be strictly increasing");
      }
    }

    if (config.mSpecies.empty()) {
      throw std::invalid_argument("Ditto::TuneAccumulator: species list must not be empty");
    }

    std::unordered_set<int> uniqueSpecies;
    uniqueSpecies.reserve(config.mSpecies.size());
    for (const int pdg : config.mSpecies) {
      if (!uniqueSpecies.insert(pdg).second) {
        throw std::invalid_argument("Ditto::TuneAccumulator: duplicate PDG code " +
                                    std::to_string(pdg));
      }
    }

    if (config.mMaxNch <= 0 || config.mMaxSelectedMultiplicity <= 0 ||
        config.mMaxSpeciesMultiplicity <= 0) {
      throw std::invalid_argument(
        "Ditto::TuneAccumulator: multiplicity histogram maxima must be > 0");
    }

    if (config.mNPtBins <= 0 || config.mPtMax <= 0.0) {
      throw std::invalid_argument("Ditto::TuneAccumulator: invalid pT binning");
    }
    if (config.mNEtaBins <= 0) {
      throw std::invalid_argument("Ditto::TuneAccumulator: nEtaBins must be > 0");
    }

    if (config.mMaxCompositionTemplatesPerPair == 0) {
      throw std::invalid_argument(
        "Ditto::TuneAccumulator: maxCompositionTemplatesPerPair must be > 0");
    }
  }

  void bookTune(const std::vector<ParticleDefinition>& particleDefinitions)
  {
    std::unordered_map<int, const ParticleDefinition*> definitionByPdg;
    definitionByPdg.reserve(particleDefinitions.size());

    for (const auto& definition : particleDefinitions) {
      if (!definitionByPdg.emplace(definition.mPdg, &definition).second) {
        throw std::invalid_argument(
          "Ditto::TuneAccumulator: duplicate particle definition for PDG code " +
          std::to_string(definition.mPdg));
      }
    }

    tune.initialize(config.mActivityEtaMax,
                    config.mParticleEtaMax,
                    config.mActivityEdges,
                    config.mSelectedMultiplicityEdges,
                    config.mMaxNch,
                    config.mMaxSelectedMultiplicity,
                    config.mMaxSpeciesMultiplicity,
                    config.mNPtBins,
                    config.mPtMax,
                    config.mNEtaBins,
                    config.mSpecies);

    tune.mAzimuthModel = "uniform";
    tune.mFinalStatus = config.mFinalStatus;

    speciesIndex.clear();
    speciesIndex.reserve(config.mSpecies.size());

    for (std::size_t i = 0; i < config.mSpecies.size(); ++i) {
      const int pdg = config.mSpecies[i];
      const auto found = definitionByPdg.find(pdg);
      if (found == definitionByPdg.end()) {
        throw std::invalid_argument(
          "Ditto::TuneAccumulator: missing particle definition for PDG code " +
          std::to_string(pdg));
      }

      speciesIndex[pdg] = i;

      const auto& definition = *found->second;
      auto* entry = tune.speciesAt(static_cast<int>(i));
      entry->mParticleName = definition.mName;
      entry->mMass = definition.mMass;
      entry->mChargeType = definition.mChargeType;
    }
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

    if (reservoir.templates.size() < config.mMaxCompositionTemplatesPerPair) {
      destination = reservoir.templates.size();
      reservoir.templates.emplace_back();
    } else {
      std::uniform_int_distribution<std::uint64_t> distribution(0, reservoir.eventsSeen - 1);
      const std::uint64_t selected = distribution(compositionRng);

      if (selected >= config.mMaxCompositionTemplatesPerPair) {
        return;
      }
      destination = static_cast<std::size_t>(selected);
    }

    auto& composition = reservoir.templates[destination];
    composition.centralCounts.resize(centralCounts.size());
    composition.otherCounts.resize(otherCounts.size());

    for (std::size_t i = 0; i < centralCounts.size(); ++i) {
      if (centralCounts[i] < 0 || otherCounts[i] < 0 ||
          centralCounts[i] > static_cast<int>(std::numeric_limits<unsigned short>::max()) ||
          otherCounts[i] > static_cast<int>(std::numeric_limits<unsigned short>::max())) {
        throw std::runtime_error(
          "Ditto::TuneAccumulator: composition multiplicity exceeds unsigned-short storage");
      }

      composition.centralCounts[i] = static_cast<unsigned short>(centralCounts[i]);
      composition.otherCounts[i] = static_cast<unsigned short>(otherCounts[i]);
    }
  }

  void serializeCompositionTemplates()
  {
    tune.mCompositionTemplateCapPerPair =
      static_cast<std::uint64_t>(config.mMaxCompositionTemplatesPerPair);

    tune.mCompositionPairNch.clear();
    tune.mCompositionPairNSelected.clear();
    tune.mCompositionPairOffsets.clear();
    tune.mCompositionPairEventsSeen.clear();
    tune.mCompositionCentralCounts.clear();
    tune.mCompositionOtherCounts.clear();

    std::vector<std::uint64_t> keys;
    keys.reserve(compositionReservoirs.size());

    for (const auto& item : compositionReservoirs) {
      if (!item.second.templates.empty()) {
        keys.push_back(item.first);
      }
    }

    std::sort(keys.begin(), keys.end());

    tune.mCompositionPairNch.reserve(keys.size());
    tune.mCompositionPairNSelected.reserve(keys.size());
    tune.mCompositionPairEventsSeen.reserve(keys.size());
    tune.mCompositionPairOffsets.reserve(keys.size() + 1);
    tune.mCompositionPairOffsets.push_back(0);

    const std::size_t nSpecies = config.mSpecies.size();
    std::uint64_t templateOffset = 0;

    for (const std::uint64_t key : keys) {
      const auto found = compositionReservoirs.find(key);
      if (found == compositionReservoirs.end()) {
        continue;
      }

      const auto& reservoir = found->second;
      tune.mCompositionPairNch.push_back(compositionKeyNch(key));
      tune.mCompositionPairNSelected.push_back(compositionKeyNSelected(key));
      tune.mCompositionPairEventsSeen.push_back(reservoir.eventsSeen);

      for (const auto& composition : reservoir.templates) {
        if (composition.centralCounts.size() != nSpecies ||
            composition.otherCounts.size() != nSpecies) {
          throw std::runtime_error(
            "Ditto::TuneAccumulator: inconsistent composition-template species count");
        }

        tune.mCompositionCentralCounts.insert(tune.mCompositionCentralCounts.end(),
                                              composition.centralCounts.begin(),
                                              composition.centralCounts.end());
        tune.mCompositionOtherCounts.insert(tune.mCompositionOtherCounts.end(),
                                            composition.otherCounts.begin(),
                                            composition.otherCounts.end());
      }

      templateOffset += static_cast<std::uint64_t>(reservoir.templates.size());
      tune.mCompositionPairOffsets.push_back(templateOffset);
    }
  }

  void processEvent(const InputEvent& event)
  {
    if (isFinalized) {
      throw std::runtime_error(
        "Ditto::TuneAccumulator: cannot process events after finalize()");
    }

    int nch = 0;

    // First pass: define event activity from all charged input particles.
    for (const auto& particle : event) {
      if (!std::isfinite(particle.mEta) ||
          std::abs(particle.mEta) >= config.mActivityEtaMax) {
        continue;
      }

      if (particle.mChargeType != 0) {
        ++nch;
      }
    }

    tune.mHNch.Fill(nch);

    if (tune.activityClass(nch) < 0) {
      ++activityOverflowEvents;
      ++processedEvents;
      return;
    }

    std::vector<int> speciesCounts(config.mSpecies.size(), 0);
    std::vector<int> centralChargedSpeciesCounts(config.mSpecies.size(), 0);
    std::vector<int> otherSpeciesCounts(config.mSpecies.size(), 0);

    int nSelected = 0;
    int nSelectedCentralCharged = 0;

    // Second pass: composition and single-particle kinematics conditioned on Nch.
    for (const auto& particle : event) {
      const auto found = speciesIndex.find(particle.mPdg);
      if (found == speciesIndex.end()) {
        continue;
      }

      const double eta = particle.mEta;
      if (!std::isfinite(eta) || std::abs(eta) >= config.mParticleEtaMax) {
        continue;
      }

      const std::size_t iSpecies = found->second;
      ++speciesCounts[iSpecies];
      ++nSelected;

      const bool centralCharged =
        particle.mChargeType != 0 && std::abs(eta) < config.mActivityEtaMax;

      if (centralCharged) {
        ++centralChargedSpeciesCounts[iSpecies];
        ++nSelectedCentralCharged;
      } else {
        ++otherSpeciesCounts[iSpecies];
      }

      const double pt = particle.mPt;
      if (pt >= config.mPtMax) {
        ++ptOverflowParticles;
      }

      auto* entry = tune.speciesAt(static_cast<int>(iSpecies));
      entry->mHPtVsActivity.Fill(nch, pt);
      entry->mHEtaVsActivity.Fill(nch, eta);

      if (centralCharged) {
        entry->mHCentralChargedPtSumVsNch.Fill(nch, pt);
        entry->mHCentralChargedPtCountVsNch.Fill(nch);
      } else {
        entry->mHOtherPtSumVsNch.Fill(nch, pt);
        entry->mHOtherPtCountVsNch.Fill(nch);
      }
    }

    tune.mHNSelected.Fill(nSelected);
    tune.mHNchSelected.Fill(nSelectedCentralCharged);

    const double selectedMultiplicityPoint[2] = {
      static_cast<double>(nch), static_cast<double>(nSelected)};
    const double selectedChargedPoint[2] = {
      static_cast<double>(nch), static_cast<double>(nSelectedCentralCharged)};

    tune.mHNSelectedVsNch.Fill(selectedMultiplicityPoint);
    tune.mHNchSelectedVsNch.Fill(selectedChargedPoint);

    const int selectedClass = tune.selectedMultiplicityClass(nSelected);
    if (selectedClass < 0) {
      ++tune.mSelectedMultiplicityOverflowEvents;
    } else {
      tune.mHEventsVsActivitySelected.Fill(nch, nSelected);
    }

    if (nSelectedCentralCharged != nch) {
      ++tune.mCentralChargedCoverageMismatchEvents;
      if (nSelectedCentralCharged < nch) {
        tune.mCentralChargedCoverageMissingParticles +=
          static_cast<std::uint64_t>(nch - nSelectedCentralCharged);
      }
    }

    // Fill once per event for every species, including N_species = 0.
    for (std::size_t i = 0; i < speciesCounts.size(); ++i) {
      const int count = speciesCounts[i];
      if (count > config.mMaxSpeciesMultiplicity) {
        ++speciesMultiplicityOverflowEvents;
      }

      auto* entry = tune.speciesAt(static_cast<int>(i));
      entry->mHCountVsActivity.Fill(nch, count);

      const int centralCount = centralChargedSpeciesCounts[i];
      const int otherCount = otherSpeciesCounts[i];
      const double centralMultiplicity[2] = {
        static_cast<double>(nch), static_cast<double>(centralCount)};
      const double otherMultiplicity[2] = {
        static_cast<double>(nch), static_cast<double>(otherCount)};

      entry->mHCentralChargedMultiplicityVsActivity.Fill(centralMultiplicity);
      entry->mHOtherMultiplicityVsActivity.Fill(otherMultiplicity);

      if (selectedClass >= 0) {
        if (centralCount > 0) {
          entry->mHCentralChargedCountVsActivitySelected.Fill(
            nch, nSelected, centralCount);
        }
        if (otherCount > 0) {
          entry->mHOtherCountVsActivitySelected.Fill(nch, nSelected, otherCount);
        }
      }
    }

    if (selectedClass >= 0) {
      storeCompositionTemplate(
        nch, nSelected, centralChargedSpeciesCounts, otherSpeciesCounts);
    }

    ++processedEvents;
  }

  void finalize()
  {
    if (isFinalized) {
      throw std::runtime_error("Ditto::TuneAccumulator: finalize() called more than once");
    }

    tune.mNEvents = processedEvents;
    tune.mActivityOverflowEvents = activityOverflowEvents;
    tune.mPtOverflowParticles = ptOverflowParticles;
    tune.mSpeciesMultiplicityOverflowEvents = speciesMultiplicityOverflowEvents;

    serializeCompositionTemplates();
    tune.finalize();
    isFinalized = true;
  }
};

TuneAccumulator::TuneAccumulator(
  const TuneAccumulatorConfig& config,
  const std::vector<ParticleDefinition>& particleDefinitions)
{
  auto impl = std::make_unique<TuneAccumulatorImpl>(config, particleDefinitions);
  mImpl = impl.release();
}

TuneAccumulator::~TuneAccumulator()
{
  delete mImpl;
  mImpl = nullptr;
}

void TuneAccumulator::processEvent(const InputEvent& event)
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::TuneAccumulator: invalid implementation");
  }
  mImpl->processEvent(event);
}

void TuneAccumulator::finalize()
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::TuneAccumulator: invalid implementation");
  }
  mImpl->finalize();
}

const TuneAccumulatorConfig& TuneAccumulator::config() const
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::TuneAccumulator: invalid implementation");
  }
  return mImpl->config;
}

std::uint64_t TuneAccumulator::processedEvents() const
{
  return mImpl ? mImpl->processedEvents : 0;
}

std::uint64_t TuneAccumulator::activityOverflowEvents() const
{
  return mImpl ? mImpl->activityOverflowEvents : 0;
}

std::uint64_t TuneAccumulator::ptOverflowParticles() const
{
  return mImpl ? mImpl->ptOverflowParticles : 0;
}

std::uint64_t TuneAccumulator::speciesMultiplicityOverflowEvents() const
{
  return mImpl ? mImpl->speciesMultiplicityOverflowEvents : 0;
}

bool TuneAccumulator::finalized() const
{
  return mImpl && mImpl->isFinalized;
}

Tune& TuneAccumulator::tune()
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::TuneAccumulator: invalid implementation");
  }
  return mImpl->tune;
}

const Tune& TuneAccumulator::tune() const
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::TuneAccumulator: invalid implementation");
  }
  return mImpl->tune;
}

} // namespace Ditto
