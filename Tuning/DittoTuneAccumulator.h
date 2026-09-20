///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTuneAccumulator.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/19
/// \brief  Generator-independent accumulation of Ditto tuning observables.
///

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Ditto
{

class Tune;

/// Generator-independent configuration of the tune accumulation.
///
/// Teacher/source-specific settings (number of generated events, PYTHIA card,
/// AO2D input, etc.) deliberately do not belong here.
struct TuneAccumulatorConfig {
  /// Event activity:
  ///   Nch = number of charged input particles in |eta| < mActivityEtaMax.
  double mActivityEtaMax = 0.5;

  /// Only particles inside this acceptance are used to build the per-species
  /// composition and kinematic tables.
  double mParticleEtaMax = 5.0;

  /// Final-state status written by Ditto when generating from the tune.
  int mFinalStatus = 1;

  std::vector<double> mActivityEdges = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0,
                                        60.0, 80.0, 100.0, 150.0, 250.0, 500.0, 1000.0};

  std::vector<double> mSelectedMultiplicityEdges = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 300.0, 500.0, 1000.0, 2000.0, 5000.0};

  int mMaxNch = 1000;
  int mMaxSelectedMultiplicity = 5000;
  int mMaxSpeciesMultiplicity = 1000;

  int mNPtBins = 1000;
  double mPtMax = 20.0;
  int mNEtaBins = 400;

  /// Exact PDG species learned separately by the tune.
  std::vector<int> mSpecies = {22,
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

  std::size_t mMaxCompositionTemplatesPerPair = 128;
  std::uint64_t mCompositionReservoirSeed = 1;
};

/// Minimal particle representation consumed by the accumulator.
///
/// The source adapter is responsible for selecting the particles that should
/// be considered final before they are passed here.
struct InputParticle {
  int mPdg = 0;
  double mPt = 0.0;
  double mEta = 0.0;
  int mChargeType = 0; ///< Three times the electric charge.
};

using InputEvent = std::vector<InputParticle>;

/// Static particle properties needed to persist a self-contained Ditto tune.
struct ParticleDefinition {
  int mPdg = 0;
  std::string mName;
  double mMass = 0.0;
  int mChargeType = 0; ///< Three times the electric charge.
};

struct TuneAccumulatorImpl;

/// Accumulates source-independent event information into a Ditto::Tune.
class TuneAccumulator
{
 public:
  TuneAccumulator(const TuneAccumulatorConfig& config,
                  const std::vector<ParticleDefinition>& particleDefinitions);
  ~TuneAccumulator();

  TuneAccumulator(const TuneAccumulator&) = delete;
  TuneAccumulator& operator=(const TuneAccumulator&) = delete;

  void processEvent(const InputEvent& event);

  /// Serialize the composition reservoirs and compute the final conditional
  /// probabilities. No more events may be processed afterwards.
  void finalize();

  const TuneAccumulatorConfig& config() const;

  std::uint64_t processedEvents() const;
  std::uint64_t activityOverflowEvents() const;
  std::uint64_t ptOverflowParticles() const;
  std::uint64_t speciesMultiplicityOverflowEvents() const;

  bool finalized() const;

  Tune& tune();
  const Tune& tune() const;

 private:
  TuneAccumulatorImpl* mImpl = nullptr;
};

} // namespace Ditto
