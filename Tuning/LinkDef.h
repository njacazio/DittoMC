///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   LinkDef.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Link definitions for Ditto.
///


#ifdef __CLING__
#pragma link C++ namespace Ditto;
#pragma link C++ class Ditto::TuneSpecies + ;
#pragma link C++ class Ditto::Tune + ;
#endif
