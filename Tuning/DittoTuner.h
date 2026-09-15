///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTuner.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Public interface for the Ditto teacher tuner.
///

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Ditto
{

class Tune;

struct TunerConfig {
  // Number of successfully generated PYTHIA events used to build the tune.
  std::uint64_t nEvents = 1000000;

  // Complete PYTHIA configuration card. All PYTHIA settings, including the
  // beams, energy, process, and random seed, must be defined in this card.
  // Ditto currently requires Beams:frameType = 1 and persists the beam IDs
  // and Beams:eCM as the generator energy metadata.
  std::string pythiaCard;

  // Event activity used for all conditional distributions:
  //
  //   Nch = number of final charged PYTHIA particles in
  //         |eta| < activityEtaMax.
  //
  double activityEtaMax = 0.5;

  // Only particles inside this acceptance are used to build the per-species
  // composition and kinematic tables.
  double particleEtaMax = 5.0;

  // Final-state status written by Ditto for particles generated from this
  // tune / generator card.
  int finalStatus = 1;

  // Nch classes used for:
  //   P(N_species | Nch class)
  //   P(pT        | species, Nch class)
  //   P(eta       | species, Nch class)
  //
  // The last edge should be comfortably above the largest expected Nch.
  std::vector<double> activityEdges = {
    0.0, 5.0, 10.0, 20.0, 30.0, 40.0,
    60.0, 80.0, 100.0, 150.0, 250.0, 500.0, 1000.0};

  // Nselected classes used only for the constrained species-composition
  // tables. The exact P(Nselected | Nch) relation is stored separately in a
  // sparse two-dimensional table.
  std::vector<double> selectedMultiplicityEdges = {
    0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 160.0,
    200.0, 300.0, 500.0, 1000.0, 2000.0, 5000.0};

  // Global histogram ranges.
  int maxNch = 1000;
  int maxSelectedMultiplicity = 5000;

  // Per-species multiplicity histogram range. This is also the exact
  // multiplicity-axis range of the sparse per-component diagnostic PMFs.
  int maxSpeciesMultiplicity = 1000;

  // Kinematic table binning.
  int nPtBins = 1000;
  double ptMax = 20.0;

  // In addition to the coarse P(pT | species, activity class), the tuner
  // stores exact-Nch mean pT separately for every species and for the
  // central-charged / other generator components. No additional binning
  // parameter is required: the exact Nch axis uses maxNch.
  int nEtaBins = 400;

  // Exact PDG species to learn separately. The list becomes part of the
  // persistent tune / generator card and is therefore also the generator list.
  std::vector<int> species;

  // Maximum number of correlated composition templates retained for each
  // exact (Nch, Nselected) pair. Reservoir sampling keeps these templates
  // uniformly distributed over all teacher events in the pair.
  std::size_t maxCompositionTemplatesPerPair = 128;
  std::uint64_t compositionReservoirSeed = 1;

  // Progress / diagnostics.
  std::uint64_t progressEvery = 10000;
  bool printPythiaStatistics = true;

  // Protect against an infinite loop if PYTHIA repeatedly fails.
  double maxAttemptsFactor = 5.0;

  TunerConfig();
};

// Opaque runtime implementation. Its definition lives only in
// DittoTuner.cxx and is intentionally not persisted.
struct TunerImpl;

/// Runs the Ditto teacher-tuning workflow for a PYTHIA configuration.
///
/// The tuner builds the tune from generated events, records the learned
/// species and composition distributions, and persists the result as a
/// `Ditto::Tune` object.
class Tuner
{
 public:
  /// Construct a tuner using the supplied configuration.
  explicit Tuner(const TunerConfig& config = TunerConfig{});

  /// Destroy the tuner and release any owned runtime state.
  ~Tuner();

  Tuner(const Tuner&) = delete;
  Tuner& operator=(const Tuner&) = delete;

  /// Run the tuning loop until the requested number of successful events has
  /// been collected or the configured maximum number of attempts is exceeded.
  void run();

  /// Return the configuration used to initialize this tuner.
  const TunerConfig& config() const;

  /// Return the number of successfully generated events used in the tune.
  std::uint64_t generatedEvents() const;

  /// Return the total number of attempted events, including failed iterations.
  std::uint64_t attemptedEvents() const;

  /// Return the learned tune produced by the tuning run.
  const Tune& tune() const;

 private:
  // Owned by Tuner. Kept as an opaque raw pointer in the public header
  // specifically to avoid ROOT dictionary generation for implementation
  // details. The constructor/destructor provide RAII ownership.
  TunerImpl* mImpl = nullptr; //! transient implementation detail
};

} // namespace Ditto
