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
/// This is deliberately a plain ROOT reader. It does not depend on DPL or on
/// O2Physics analysis tasks. Input strings are passed to TFile::Open(), so
/// ROOT-supported remote protocols (including alien:// when available in the
/// running environment) work transparently.
struct AO2DTunerConfig : public TuneAccumulatorConfig {
  /// Input AO2D files. Local paths and ROOT-supported remote URLs are accepted.
  std::vector<std::string> inputFiles;

  /// MC-particle selection used to define the source event.
  AO2DParticleSelection particleSelection = AO2DParticleSelection::GeneratorFinal;

  /// Stop after this many MC collisions. 0 means process every input event.
  std::uint64_t maxEvents = 0;

  /// Print a progress line every N processed events. 0 disables progress.
  std::uint64_t progressEvery = 10000;

  /// Beam metadata to persist in the tune. AO2D itself does not provide these
  /// in a source-independent way, so they are explicit configuration here.
  int beamIdA = 0;
  int beamIdB = 0;
  int beamFrameType = 1;
  double sqrtSNN = 0.0;

  /// Unknown PDG codes make Nch ambiguous because the charge is unknown.
  /// Therefore the default is to fail loudly. If enabled, unknown particles
  /// are skipped entirely and a warning count is printed at the end.
  bool ignoreUnknownPdg = false;
};

struct AO2DTunerImpl;

/// Standalone AO2D source adapter for TuneAccumulator.
///
/// Typical ROOT-macro usage:
///
///   Ditto::AO2DTunerConfig cfg;
///   cfg.inputFiles = {"AO2D.root"};
///   cfg.beamIdA = 2212;
///   cfg.beamIdB = 2212;
///   cfg.sqrtSNN = 13600.;
///
///   Ditto::AO2DTuner tuner(cfg);
///   tuner.run();
///   tuner.save("Ditto_tune_AO2D.root");
///
class AO2DTuner
{
 public:
  explicit AO2DTuner(const AO2DTunerConfig& config = AO2DTunerConfig{});
  ~AO2DTuner();

  AO2DTuner(const AO2DTuner&) = delete;
  AO2DTuner& operator=(const AO2DTuner&) = delete;

  /// Add an input after construction but before run().
  void addFile(const std::string& fileName);

  /// Process all configured input files (or maxEvents) and finalize the tune.
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
