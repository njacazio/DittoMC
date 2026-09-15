///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   rootlogon.C
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  ROOT startup configuration for Ditto.
///

{

  TString dittomcRoot = gSystem->Getenv("DITTOMC_ROOT");
  if (dittomcRoot.IsNull()) {
    gSystem->AddIncludePath("-I./include/Ditto/");
    gSystem->AddDynamicPath("./lib/");
    gSystem->Load("./lib/libDitto.so");
    gSystem->Load("./lib/libDittoTuner.so");
    Printf("Ready!");
  }
}
