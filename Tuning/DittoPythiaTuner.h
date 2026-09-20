///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoPythiaTuner.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Public interface for the Ditto PYTHIA teacher tuner.
///

#pragma once

#include "DittoTuneAccumulator.h"

#include <cstdint>
#include <string>

namespace Ditto
{

class Tune;

/// PYTHIA-specific tuning configuration.
///
/// The generator-independent accumulation settings are inherited from
/// TuneAccumulatorConfig so existing code such as cfg.mActivityEtaMax remains
/// unchanged.
struct TunerConfig : public TuneAccumulatorConfig {
  /// Number of successfully generated PYTHIA events used to build the tune.
  std::uint64_t nEvents = 1000000;

  /// Complete PYTHIA configuration card. All PYTHIA settings, including the
  /// beams, energy, process, and random seed, must be defined in this card.
  /// Ditto currently requires Beams:frameType = 1 and persists the beam IDs
  /// and Beams:eCM as the generator energy metadata.
  std::string pythiaCard;

  /// Progress / diagnostics.
  std::uint64_t progressEvery = 10000;
  bool printPythiaStatistics = true;

  /// Protect against an infinite loop if PYTHIA repeatedly fails.
  double maxAttemptsFactor = 5.0;
};

// Opaque runtime implementation. Its definition lives only in
// DittoPythiaTuner.cxx and is intentionally not persisted.
struct TunerImpl;

/// PYTHIA source adapter for TuneAccumulator.
class PythiaTuner
{
 public:
  /// Construct a tuner using the supplied configuration.
  explicit PythiaTuner(const TunerConfig& config = TunerConfig{});

  /// Destroy the tuner and release any owned runtime state.
  ~PythiaTuner();

  PythiaTuner(const PythiaTuner&) = delete;
  PythiaTuner& operator=(const PythiaTuner&) = delete;

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
  // Owned by PythiaTuner. Kept as an opaque raw pointer in the public header
  // specifically to avoid ROOT dictionary generation for implementation
  // details. The constructor/destructor provide RAII ownership.
  TunerImpl* mImpl = nullptr; //! transient implementation detail
};

} // namespace Ditto
