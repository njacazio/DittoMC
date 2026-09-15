///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   exampleTuner.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Example Ditto tuning macro.
///

#include "DittoTuner.h"

#include <cstdint>
#include <iostream>

void exampleTuner(std::uint64_t nEvents, const std::string pythiaCard)
{
  Ditto::TunerConfig cfg;

  cfg.pythiaCard = pythiaCard.c_str();
  cfg.nEvents = nEvents;

  // Event activity for the conditional tables.
  cfg.activityEtaMax = 0.5;

  // Ditto generation acceptance that we want to learn.
  cfg.particleEtaMax = 5.0;

  // This also becomes part of the persistent generator card.
  cfg.finalStatus = 1;

  if (pythiaCard.find("pythia8_inel_136tev") != std::string::npos) { // pp

    // Enough for pp at 13 TeV. Increase these for PbPb.
    cfg.activityEdges = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0,
                         60.0, 80.0, 100.0, 150.0};

    // Coarser classes are sufficient for the species-composition fractions.
    // P(Nselected | Nch) itself is stored with exact integer multiplicities.
    cfg.selectedMultiplicityEdges = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0,
                                     120.0, 160.0, 200.0, 300.0, 500.0, 750.0,
                                     1000.0};
    cfg.maxNch = 200;
    cfg.maxSelectedMultiplicity = 1000;
    // Also controls the z range of the diagnostic species-multiplicity PMFs.
    cfg.maxSpeciesMultiplicity = 500;
    cfg.progressEvery = 10000;
  } else if (pythiaCard.find("pythia8_OO_536") != std::string::npos) { // OO
    cfg.activityEdges = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0,
                         60.0, 80.0, 100.0, 150.0, 250.0, 500.0, 1000.0};

    // Coarser classes are sufficient for the species-composition fractions.
    // P(Nselected | Nch) itself is stored with exact integer multiplicities.
    cfg.selectedMultiplicityEdges = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0,
                                     120.0, 160.0, 200.0, 300.0, 500.0, 750.0,
                                     1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0, 10000.0};
    cfg.maxNch = 2000;
    cfg.maxSelectedMultiplicity = 5000;
    // Also controls the z range of the diagnostic species-multiplicity PMFs.
    cfg.maxSpeciesMultiplicity = 1000;
    cfg.progressEvery = 10;
  } else if (pythiaCard.find("pythia8_PbPb_536tev") != std::string::npos) { // PbPb

    cfg.activityEdges = {0., 50., 100., 200., 300., 500.,
                         750., 1000., 1500., 2000., 3000.,
                         4000., 5000.};

    cfg.selectedMultiplicityEdges = {0., 100., 200., 500., 1000., 2000.,
                                     5000., 10000., 20000., 30000., 50000.};

    cfg.maxNch = 5000;
    cfg.maxSelectedMultiplicity = 50000;
    cfg.maxSpeciesMultiplicity = 10000;

    cfg.progressEvery = 10;
  } else {
    throw std::invalid_argument("exampleTuner: unrecognized PYTHIA card name: " + pythiaCard);
  }

  // Uniform reservoir of full correlated composition vectors for each exact
  // (Nch, Nselected) pair.
  cfg.maxCompositionTemplatesPerPair = 128;
  cfg.compositionReservoirSeed = 1;

  cfg.ptMax = 20.0;
  cfg.nPtBins = 1000;
  cfg.nEtaBins = 400;

  Ditto::Tuner tuner(cfg);
  tuner.run();

  std::cout << "Generated " << tuner.generatedEvents() << " successful PYTHIA events\n";
}

void exampleTuner(int cfg = 0)
{

  std::uint64_t nEvents = 0;
  std::string pythiaCard = "";
  switch (cfg) {
    case 0:
      nEvents = 1E5;
      pythiaCard = "Tuning/cards/pythia8_inel_136tev.cfg";
      break;
    case 1:
      nEvents = 1E6;
      pythiaCard = "Tuning/cards/pythia8_OO_536.cfg";
      break;
    case 2:
      nEvents = 1E4;
      pythiaCard = "Tuning/cards/pythia8_PbPb_536tev.cfg";
      break;
    default:
      throw std::invalid_argument("exampleTuner: unrecognized cfg value: " + std::to_string(cfg));
  }
  exampleTuner(nEvents, pythiaCard);
}
