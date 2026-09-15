///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   Ditto.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Public interface for the Ditto generator.
///

#pragma once

#include <TClonesArray.h>
#include <TParticle.h>

#include <Pythia8/Pythia.h>

#include <cstdint>
#include <iosfwd>
#include <random>
#include <string>
#include <string_view>
#include <vector>

class TFile;
class TTree;

namespace Ditto
{

// Runtime options only.
//
// All physics configuration belongs to the tune file: beam IDs, energy,
// generated acceptance, particle species and the learned event/kinematic model.
struct Config {
  // Ditto tune / generator-card file produced by Ditto::Tuner.
  std::string tuneFile;

  // Ditto random seed. This is intentionally independent of the PYTHIA
  // random seed stored in the teacher card.
  std::uint64_t seed = 1;

  // Measure total generator timing and print a summary when the Generator is
  // destroyed. Pure generation, TTree filling and PYTHIA export are reported
  // separately.
  bool enableTimingMetrics = false;

  // Add fine-grained timers inside generation. This is intended for profiling:
  // it performs many steady_clock::now() calls per particle and therefore
  // perturbs the absolute timing. Use enableTimingMetrics=true and this=false
  // for the least intrusive throughput benchmark.
  bool enableDetailedTimingMetrics = false;
};

struct EventInfo {
  std::uint64_t eventNumber = 0;

  // Teacher beam configuration carried by the tune / generator card.
  int beamIdA = 0;
  int beamIdB = 0;
  double sqrtSNN = 0.0;

  // Exact event-level conditioning variables sampled from the tune.
  int conditioningNch = -1;
  int conditioningNSelected = -1;
  int activityClass = -1;
};

// Accumulated wall-clock timings. Values are stored in seconds.
//
// The detailed particle components are only populated when
// Config::enableDetailedTimingMetrics is true.
struct TimingMetrics {
  std::uint64_t generatedEvents = 0;
  std::uint64_t generatedParticles = 0;
  std::uint64_t treeFills = 0;
  std::uint64_t pythiaExports = 0;

  // Coarse timings.
  double generation = 0.0; // excludes TTree::Fill and PYTHIA export
  double treeFill = 0.0;
  double pythiaExport = 0.0;

  // Event-level generation components.
  double clearParticles = 0.0;
  double eventInfo = 0.0;
  double multiplicity = 0.0;
  double expandArray = 0.0;
  double particleLoop = 0.0;

  // Particle-generation components.
  double speciesSampling = 0.0;
  double constructedAt = 0.0;
  double ptSampling = 0.0;
  double etaSampling = 0.0;
  double phiSampling = 0.0;
  double momentumMath = 0.0;
  double particleSetters = 0.0;

  // PYTHIA export components.
  double pythiaReset = 0.0;
  double pythiaAppend = 0.0;
  double pythiaSystemSum = 0.0;
};

class Tune;
struct TuneAliasSampler;
struct TuneDiscreteAliasSampler;
struct TuneRuntimeData;

class Generator
{
 public:
  explicit Generator(const Config& config);

  // Destructor, writes the TTree output if registered.
  ~Generator();

  // Generate and store one Ditto event internally.
  const EventInfo& generate();

  // Alias convenient in an event loop.
  const EventInfo& next() { return generate(); }

  const EventInfo& currentEvent() const { return mEvent; }
  const Config& config() const { return mConfig; }

  // Persistent tune object used as the complete physics generator card.
  const Tune& generatorCard() const;

  // Average pure generation time, excluding TTree filling and PYTHIA export.
  // Returns 0 if timing metrics are disabled or no event has been generated.
  double averageGenerationTimeUs() const;

  const TimingMetrics& timingMetrics() const { return mTiming; }
  void resetTimingMetrics();
  void printTimingMetrics(std::ostream& os) const;

  // Export the already-generated event. No new random numbers are consumed.
  // If reset=true, the target PYTHIA event is reset before filling.
  void loadParticles(Pythia8::Event& event, bool reset = true) const;

  /// @brief Creates an output file with a TTree containing the generated events
  /// @param fileName The name of the output file
  void registerTTreeOutput(const std::string& fileName);

 private:
  Config mConfig;
  EventInfo mEvent;
  std::mt19937_64 mRng;
  std::uint64_t mEventCounter = 0;

  mutable TimingMetrics mTiming;

  // The tune is mandatory and is the complete physics configuration.
  Tune* mTune = nullptr;                   //!
  TuneRuntimeData* mTuneRuntime = nullptr; //!

  void validateConfig() const;
  void loadTune(const std::string& fileName);

  bool timingEnabled() const
  {
    return mConfig.enableTimingMetrics || mConfig.enableDetailedTimingMetrics;
  }

  bool detailedTimingEnabled() const
  {
    return mConfig.enableDetailedTimingMetrics;
  }

  double uniform01();
  double uniform(double min, double max);

  EventInfo makeEventInfo() const;
  int sampleMultiplicity(EventInfo& info);

  void sampleComposition(const EventInfo& info,
                         std::vector<int>& centralCounts,
                         std::vector<int>& otherCounts);

  double sampleAlias(const TuneAliasSampler& sampler);
  int sampleDiscreteAlias(const TuneDiscreteAliasSampler& sampler);
  double samplePhi();

  void makeParticle(TParticle* particle,
                    int pdg,
                    double mass,
                    double mass2,
                    double ptScale,
                    const TuneAliasSampler& ptSampler,
                    const TuneAliasSampler& etaSampler);

  // Output file for TTree writing. If empty, no output is written.
  static constexpr std::string_view kTreeName = "T";
  static constexpr std::string_view kBranchName = "Particles";
  TFile* mOutputFile = nullptr;
  TTree* mTree = nullptr;
  TClonesArray* mParticles = nullptr;
};

} // namespace Ditto
