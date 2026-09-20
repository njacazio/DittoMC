///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoPythiaTuner.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  PYTHIA source adapter for the Ditto tune accumulator.
///

#include "DittoPythiaTuner.h"

#include "DittoTune.h"

// ACLiC/rootcling only needs to see the public declarations from the header.
// Hide all implementation details from dictionary generation.
#ifndef __ROOTCLING__

#include <Pythia8/Pythia.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
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
    throw std::runtime_error("Ditto::PythiaTuner: could not read PYTHIA card: " + fileName);
  }

  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

} // namespace

struct TunerImpl {
  explicit TunerImpl(const TunerConfig& cfg)
    : config(cfg)
  {
  }

  TunerConfig config;
  std::string outputFile;

  std::unique_ptr<Pythia8::Pythia> pythia;
  std::unique_ptr<TuneAccumulator> accumulator;
  InputEvent eventBuffer;

  std::uint64_t attemptedEvents = 0;
  bool ran = false;

  void validateSourceConfig() const
  {
    if (config.pythiaCard.empty()) {
      throw std::invalid_argument("Ditto::PythiaTuner: pythiaCard must not be empty");
    }
    if (config.nEvents == 0) {
      throw std::invalid_argument("Ditto::PythiaTuner: nEvents must be > 0");
    }
    if (config.maxAttemptsFactor < 1.0) {
      throw std::invalid_argument("Ditto::PythiaTuner: maxAttemptsFactor must be >= 1");
    }
  }

  void configurePythia()
  {
    if (!pythia->readFile(config.pythiaCard)) {
      throw std::runtime_error("Ditto::PythiaTuner: could not read PYTHIA card: " +
                               config.pythiaCard);
    }

    const int frameType = pythia->settings.mode("Beams:frameType");
    if (frameType != 1) {
      throw std::invalid_argument(
        "Ditto::PythiaTuner: only Beams:frameType = 1 is currently supported");
    }

    const double sqrtSNN = pythia->settings.parm("Beams:eCM");
    if (!std::isfinite(sqrtSNN) || sqrtSNN <= 0.0) {
      throw std::invalid_argument("Ditto::PythiaTuner: Beams:eCM must be > 0");
    }

    if (!pythia->init()) {
      throw std::runtime_error("Ditto::PythiaTuner: PYTHIA initialization failed");
    }

    for (const int pdg : config.mSpecies) {
      if (!pythia->particleData.isParticle(pdg)) {
        throw std::invalid_argument("Ditto::PythiaTuner: unknown PDG code " +
                                    std::to_string(pdg));
      }
    }
  }

  void configureAccumulator()
  {
    std::vector<ParticleDefinition> particleDefinitions;
    particleDefinitions.reserve(config.mSpecies.size());

    for (const int pdg : config.mSpecies) {
      particleDefinitions.push_back({pdg,
                                     pythia->particleData.name(pdg),
                                     pythia->particleData.m0(pdg),
                                     pythia->particleData.chargeType(pdg)});
    }

    accumulator = std::make_unique<TuneAccumulator>(config, particleDefinitions);

    auto& tune = accumulator->tune();
    tune.mTeacher = "PYTHIA8";
    tune.mPythiaCard = config.pythiaCard;
    tune.mPythiaCardContent = readTextFile(config.pythiaCard);
    tune.mBeamIdA = pythia->settings.mode("Beams:idA");
    tune.mBeamIdB = pythia->settings.mode("Beams:idB");
    tune.mBeamFrameType = pythia->settings.mode("Beams:frameType");
    tune.mSqrtSNN = pythia->settings.parm("Beams:eCM");
  }

  std::uint64_t generatedEvents() const
  {
    return accumulator ? accumulator->processedEvents() : 0;
  }

  Tune& tune()
  {
    if (!accumulator) {
      throw std::runtime_error("Ditto::PythiaTuner: accumulator is not initialized");
    }
    return accumulator->tune();
  }

  const Tune& tune() const
  {
    if (!accumulator) {
      throw std::runtime_error("Ditto::PythiaTuner: accumulator is not initialized");
    }
    return accumulator->tune();
  }

  void analyzeCurrentEvent()
  {
    const auto& event = pythia->event;

    eventBuffer.clear();
    eventBuffer.reserve(static_cast<std::size_t>(event.size()));

    for (int i = 1; i < event.size(); ++i) {
      const auto& particle = event[i];
      if (!particle.isFinal()) {
        continue;
      }

      eventBuffer.push_back({particle.id(),
                             particle.pT(),
                             particle.eta(),
                             particle.chargeType()});
    }

    accumulator->processEvent(eventBuffer);
  }

  void saveTune()
  {
    auto& outputTune = tune();
    outputTune.mNAttempts = attemptedEvents;

    accumulator->finalize();
    outputTune.save(outputFile);
  }

  void run()
  {
    if (ran) {
      throw std::runtime_error(
        "Ditto::PythiaTuner: run() may only be called once per PythiaTuner instance");
    }
    ran = true;

    const std::uint64_t maxAttempts = static_cast<std::uint64_t>(
      std::ceil(config.maxAttemptsFactor * static_cast<double>(config.nEvents)));

    const auto start = std::chrono::steady_clock::now();

    while (generatedEvents() < config.nEvents) {
      if (attemptedEvents >= maxAttempts) {
        throw std::runtime_error("Ditto::PythiaTuner: too many failed PYTHIA attempts");
      }

      ++attemptedEvents;
      if (!pythia->next()) {
        continue;
      }

      analyzeCurrentEvent();

      const auto nGenerated = generatedEvents();
      if (config.progressEvery > 0 && nGenerated % config.progressEvery == 0) {
        const double elapsed = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
        const double eventsPerSecond =
          elapsed > 0.0 ? static_cast<double>(nGenerated) / elapsed : 0.0;
        const double etaSeconds =
          eventsPerSecond > 0.0
            ? static_cast<double>(config.nEvents - nGenerated) / eventsPerSecond
            : 0.0;

        std::cout << "Ditto tuner: "
                  << nGenerated
                  << " / "
                  << config.nEvents
                  << " successful PYTHIA events (ETA: "
                  << etaSeconds
                  << " s)\n";
      }
    }

    const auto end = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(end - start).count();
    tune().mGenerationTimeSeconds = elapsed;

    const auto nGenerated = generatedEvents();

    std::cout << "\nDitto tuner generation complete\n"
              << "  PYTHIA card      : " << config.pythiaCard << "\n"
              << "  beams            : " << tune().mBeamIdA << " + " << tune().mBeamIdB << "\n"
              << "  sqrt(s_NN)       : " << tune().mSqrtSNN << " GeV\n"
              << "  events           : " << nGenerated << "\n"
              << "  attempts         : " << attemptedEvents << "\n"
              << "  wall time        : " << std::fixed << std::setprecision(3)
              << elapsed << " s\n";

    if (elapsed > 0.0) {
      std::cout << "  teacher rate     : "
                << std::setprecision(2)
                << static_cast<double>(nGenerated) / elapsed
                << " events/s\n";
    }

    if (accumulator->activityOverflowEvents() > 0) {
      std::cout << "  WARNING activity overflow events: "
                << accumulator->activityOverflowEvents()
                << " (extend activityEdges)\n";
    }

    if (tune().mSelectedMultiplicityOverflowEvents > 0) {
      std::cout << "  WARNING selected-multiplicity overflow events: "
                << tune().mSelectedMultiplicityOverflowEvents
                << " (extend selectedMultiplicityEdges)\n";
    }

    if (accumulator->ptOverflowParticles() > 0) {
      std::cout << "  WARNING pT overflow particles: "
                << accumulator->ptOverflowParticles()
                << " (increase ptMax if relevant)\n";
    }

    if (tune().mCentralChargedCoverageMismatchEvents > 0) {
      const double eventFraction =
        static_cast<double>(tune().mCentralChargedCoverageMismatchEvents) /
        static_cast<double>(nGenerated);

      std::cout << "  WARNING selected species do not cover all central charged particles in "
                << tune().mCentralChargedCoverageMismatchEvents
                << " events ("
                << std::setprecision(4)
                << 100.0 * eventFraction
                << "%), missing "
                << tune().mCentralChargedCoverageMissingParticles
                << " charged particles in total\n";
    }

    if (accumulator->speciesMultiplicityOverflowEvents() > 0) {
      std::cout << "  WARNING species-count overflow fills: "
                << accumulator->speciesMultiplicityOverflowEvents()
                << " (increase maxSpeciesMultiplicity)\n";
    }

    if (config.printPythiaStatistics) {
      pythia->stat();
    }

    std::cout << "\nComputing conditional probabilities and writing tune...\n";

    saveTune();

    std::cout << "  composition pairs     : "
              << tune().numberOfCompositionPairs()
              << "\n"
              << "  templates stored      : "
              << tune().numberOfCompositionTemplates()
              << "\n"
              << "  template cap / pair   : "
              << tune().mCompositionTemplateCapPerPair
              << "\n"
              << "Tune written to: "
              << outputFile
              << "\n";
  }
};

PythiaTuner::PythiaTuner(const TunerConfig& config)
{
  auto impl = std::make_unique<TunerImpl>(config);

  impl->validateSourceConfig();

  const std::string cardDir = cardPath(impl->config.pythiaCard);
  const std::string cardFileStem = cardStem(impl->config.pythiaCard);
  const std::string outputDir =
    (cardDir.empty() ? "." : cardDir) + "/../tunes/";

  if (!std::filesystem::exists(outputDir)) {
    std::filesystem::create_directories(outputDir);
  }

  impl->outputFile = outputDir + "Ditto_tune_" + cardFileStem + ".root";

  // Use PYTHIA's standard construction / installation lookup. In particular,
  // PYTHIA8DATA is honored by PYTHIA itself when it is set.
  impl->pythia = std::make_unique<Pythia8::Pythia>();
  impl->configurePythia();
  impl->configureAccumulator();

  mImpl = impl.release();
}

PythiaTuner::~PythiaTuner()
{
  delete mImpl;
  mImpl = nullptr;
}

void PythiaTuner::run()
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::PythiaTuner: invalid implementation");
  }
  mImpl->run();
}

const TunerConfig& PythiaTuner::config() const
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::PythiaTuner: invalid implementation");
  }
  return mImpl->config;
}

std::uint64_t PythiaTuner::generatedEvents() const
{
  return mImpl ? mImpl->generatedEvents() : 0;
}

std::uint64_t PythiaTuner::attemptedEvents() const
{
  return mImpl ? mImpl->attemptedEvents : 0;
}

const Tune& PythiaTuner::tune() const
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::PythiaTuner: invalid implementation");
  }
  return mImpl->tune();
}

} // namespace Ditto

#endif // __ROOTCLING__
