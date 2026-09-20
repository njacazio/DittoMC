///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoAO2DTuner.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/20
/// \brief  Standalone ROOT/AO2D source adapter for Ditto tune accumulation.
///

#pragma once

#include "DittoTuneAccumulator.h"

#include <cstdint>
#include <string>
#include <vector>

namespace Ditto
{

class Tune;

/// Selection applied to the stored AO2D MC particles before handing an event
/// to TuneAccumulator.
enum class AO2DParticleSelection {
  /// Closest source-independent equivalent of Pythia8::Particle::isFinal():
  /// particle is generator-produced and has HepMC status code 1.
  GeneratorFinal,

  /// ALICE physical-primary definition encoded in the AO2D MC-particle flags.
  PhysicalPrimary
};

/// Configuration for tuning directly from one or more AO2D files.
///
/// The implementation is deliberately independent of DPL and O2Physics
/// analysis tasks. The first input file is inspected to discover the exact
/// versioned AO2D table names. All DF_* particle tables are then concatenated
/// into a ROOT TChain, while DF boundaries are retained internally so that
/// local AO2D collision indices remain correct.
///
/// Local paths and ROOT-supported remote URLs are accepted. In particular,
/// alien:// inputs work when the ROOT/AliEn plugin is available.
struct AO2DTunerConfig : public TuneAccumulatorConfig {
  /// Input AO2D files. Local paths and ROOT-supported remote URLs are accepted.
  std::vector<std::string> mInputFiles;

  /// MC-particle selection used to define the source event.
  AO2DParticleSelection mParticleSelection = AO2DParticleSelection::GeneratorFinal;

  /// Stop after this many MC collisions. 0 means process every input event.
  std::uint64_t mMaxEvents = 0;

  /// Print a progress line every N processed events. 0 disables progress.
  std::uint64_t mProgressEvery = 10000;

  /// Reference PYTHIA card metadata persisted in the tune file.
  std::string mPythiaCard;

  /// Number of AO2D files processed in one TChain.
  /// 0 means all files in a single batch.
  std::size_t mFileBatchSize = 10;

  /// Unknown PDG codes make Nch ambiguous because the charge is unknown.
  /// Therefore the default is to fail loudly. If enabled, unknown particles
  /// are skipped entirely and a warning count is printed at the end.
  bool mIgnoreUnknownPdg = false;
};

struct AO2DTunerImpl;

/// Standalone AO2D source adapter for TuneAccumulator.
///
/// Typical ROOT-macro usage:
///
/// \code{.cpp}
/// Ditto::AO2DTunerConfig cfg;
/// cfg.mInputFiles = {"AO2D.root"};
/// cfg.mPythiaCard = "Tuning/cards/pythia8_inel_136tev.cfg";
///
/// Ditto::AO2DTuner tuner(cfg);
/// tuner.run();
/// tuner.save("Ditto_tune_AO2D.root");
/// \endcode
class AO2DTuner
{
 public:
  explicit AO2DTuner(const AO2DTunerConfig& config = AO2DTunerConfig{});
  ~AO2DTuner();

  AO2DTuner(const AO2DTuner&) = delete;
  AO2DTuner& operator=(const AO2DTuner&) = delete;

  /// Process all configured input files (or mMaxEvents) and finalize the tune.
  /// May only be called once per object.
  void run();

  /// Save the finalized tune. run() must have completed successfully.
  void save(const std::string& fileName) const;

  const AO2DTunerConfig& config() const;
  std::uint64_t processedEvents() const;
  std::uint64_t processedParticles() const;
  std::uint64_t selectedParticles() const;
  std::uint64_t unknownPdgParticles() const;
  const Tune& tune() const;

 private:
  AO2DTunerImpl* mImpl = nullptr; //! transient implementation detail
};

} // namespace Ditto
