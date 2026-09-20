R__LOAD_LIBRARY(libDitto.so)

#include "DittoAO2DTuner.h"

void tuneAO2D(const std::string& inputFile = "aodlist.txt")
{
  Ditto::AO2DTunerConfig cfg;

  // Read the input file list.
  std::ifstream file(inputFile);
  std::string line;
  while (std::getline(file, line)) {
    if (!line.empty()) {
      cfg.inputFiles.push_back(line);
    }
  }

  cfg.beamIdA = 2212;
  cfg.beamIdB = 2212;
  cfg.beamFrameType = 1;
  cfg.sqrtSNN = 13600.;

  // 0 = all MC collisions in all input files.
  cfg.maxEvents = 0;
  cfg.progressEvery = 10000;

  // Closest AO2D equivalent to PYTHIA isFinal().
  cfg.particleSelection = Ditto::AO2DParticleSelection::GeneratorFinal;

  // The normal TuneAccumulator configuration is available here too.
  cfg.activityEtaMax = 0.5;
  cfg.particleEtaMax = 5.0;
  cfg.ptMax = 20.0;

  Ditto::AO2DTuner tuner(cfg);
  tuner.run();
  tuner.save("Ditto_tune_AO2D.root");
}
