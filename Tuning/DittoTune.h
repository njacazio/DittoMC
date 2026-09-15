///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTune.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Persistent data model for Ditto generator tunes.
///

#ifndef Ditto_TUNE_H
#define Ditto_TUNE_H

#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TObjArray.h>
#include <TObject.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace Ditto
{

class TuneSpecies : public TObject
{
 public:
  int pdg = 0;
  std::string particleName;

  // Particle properties copied from the initialized PYTHIA teacher so the
  // Ditto generator does not need a runtime particle-data lookup.
  double mass = 0.0;
  int chargeType = 0; // PYTHIA convention: three times the electric charge.

  TH2D hCountVsActivity;
  TH2D hPtVsActivity;
  TH2D hEtaVsActivity;

  TH2D pCountGivenActivity;
  TH2D pPtGivenActivity;
  TH2D pEtaGivenActivity;

  // Exact-Nch, species-dependent pT moments. These retain the multiplicity
  // dependence washed out by the coarser activity classes used by
  // pPtGivenActivity.
  //
  // "Central charged" and "other" use the same component definition as the
  // empirical composition templates. The pT-sum histograms have Sumw2 enabled:
  // their squared bin errors therefore store sum(pT^2), allowing the finalized
  // mean-pT histograms to carry the statistical error on the mean.
  TH1D hCentralChargedPtSumVsNch;
  TH1D hCentralChargedPtCountVsNch;
  TH1D hOtherPtSumVsNch;
  TH1D hOtherPtCountVsNch;

  TH1D hCentralChargedMeanPtVsNch;
  TH1D hOtherMeanPtVsNch;

  // Particle-count accumulators used to learn a constrained event
  // composition. "Central charged" means a charged selected particle inside
  // |eta| < activityEtaMax. "Other" contains every other selected particle.
  //
  // The axes are:
  //   x = Nch activity class
  //   y = Nselected class
  TH2D hCentralChargedCountVsActivitySelected;
  TH2D hOtherCountVsActivitySelected;

  // Species probabilities normalized across all species in each
  // (activity, Nselected) cell. Retained as useful composition diagnostics.
  TH2D pCentralChargedSpeciesGivenActivitySelected;
  TH2D pOtherSpeciesGivenActivitySelected;

  // Full per-species multiplicity distributions retained as diagnostics.
  // Generation uses the empirical correlated composition templates below.
  // These are filled once per event, including zero multiplicity. THnSparseD
  // keeps the large exact-multiplicity axis compact for heavy-ion tunes.
  //
  // Axes:
  //   0 = Nch activity class
  //   1 = exact species multiplicity
  THnSparseD hCentralChargedMultiplicityVsActivity;
  THnSparseD hOtherMultiplicityVsActivity;

  THnSparseD pCentralChargedMultiplicityGivenActivity;
  THnSparseD pOtherMultiplicityGivenActivity;

  TuneSpecies();
  TuneSpecies(int pdg,
              int nActivityBins,
              const double* activityEdges,
              int nSelectedBins,
              const double* selectedMultiplicityEdges,
              int maxSpeciesMultiplicity,
              int maxNch,
              int nPtBins,
              double ptMax,
              int nEtaBins,
              double etaMax);

  void finalize();
  void detachFromDirectories();

  ClassDefOverride(TuneSpecies, 6);
};

// Persistent Ditto generator card.
//
// It contains both teacher provenance and every physics input required by the
// runtime generator. A Generator never accepts an independent beam energy,
// collision system, acceptance or particle list.
class Tune : public TObject
{
 public:
  static constexpr const char* kObjectName = "DittoTune";
  static constexpr int kFormatVersion = 8;

  int formatVersion = kFormatVersion;
  bool finalized = false;

  std::string teacher = "PYTHIA8";
  std::string pythiaCard;
  std::string pythiaCardContent;

  // Complete beam / generator metadata copied from the initialized teacher.
  // Ditto does not accept an independent collision system or energy: these
  // values are the generator configuration.
  int beamIdA = 0;
  int beamIdB = 0;
  int beamFrameType = 0;
  double sqrtSNN = 0.0;

  // Current Ditto kinematic model samples the absolute azimuth uniformly.
  // Keeping this explicit in the card prevents an implicit analytic fallback.
  std::string azimuthModel = "uniform";

  // Final-state status written to TParticle / exported PYTHIA events.
  int finalStatus = 1;

  std::uint64_t nEvents = 0;
  std::uint64_t nAttempts = 0;
  double generationTimeSeconds = 0.0;

  std::uint64_t activityOverflowEvents = 0;
  std::uint64_t selectedMultiplicityOverflowEvents = 0;
  std::uint64_t ptOverflowParticles = 0;
  std::uint64_t speciesMultiplicityOverflowEvents = 0;

  // Diagnostic for the exact Nch constraint. These count events/particles
  // for which final charged PYTHIA particles in the activity acceptance are
  // not represented by the selected species list.
  std::uint64_t centralChargedCoverageMismatchEvents = 0;
  std::uint64_t centralChargedCoverageMissingParticles = 0;

  double activityEtaMax = 0.0;
  double particleEtaMax = 0.0;
  double ptMax = 0.0;

  std::vector<double> activityEdges;
  std::vector<double> selectedMultiplicityEdges;

  TH1D hNch;
  TH1D hNSelected;
  TH1D hNchSelected;

  // Exact event-level correlations. PbPb occupies only a tiny fraction of the
  // full integer (Nch, Nselected) plane, so these are sparse. The raw
  // accumulators are tuner-only; the normalized conditional distributions
  // below are persisted and used by the generator.
  THnSparseD hNSelectedVsNch; //!
  THnSparseD hNchSelectedVsNch; //!

  // Event occupancy for the coarser composition conditioning.
  TH2D hEventsVsActivitySelected;

  TH1D pNch;
  TH1D pNSelected;
  TH1D pNchSelected;

  THnSparseD pNSelectedGivenNch;
  THnSparseD pNchSelectedGivenNch;

  // Empirical correlated species-composition templates.
  //
  // Pairs are stored in lexicographic (Nch, Nselected) order. For pair i,
  // templates are in the half-open range
  //
  //   [compositionPairOffsets[i], compositionPairOffsets[i + 1])
  //
  // and each template contains numberOfSpecies() central and other counts in
  // the two flattened arrays below.
  std::uint64_t compositionTemplateCapPerPair = 0;

  std::vector<int> compositionPairNch;
  std::vector<int> compositionPairNSelected;
  std::vector<std::uint64_t> compositionPairOffsets;
  std::vector<std::uint64_t> compositionPairEventsSeen;

  std::vector<unsigned short> compositionCentralCounts;
  std::vector<unsigned short> compositionOtherCounts;

  // Owned TuneSpecies objects.
  TObjArray species;

  Tune();
  Tune(const Tune& other);
  Tune& operator=(const Tune& other);
  ~Tune() override;

  void initialize(double activityEtaMax,
                  double particleEtaMax,
                  const std::vector<double>& activityEdges,
                  const std::vector<double>& selectedMultiplicityEdges,
                  int maxNch,
                  int maxSelectedMultiplicity,
                  int maxSpeciesMultiplicity,
                  int nPtBins,
                  double ptMax,
                  int nEtaBins,
                  const std::vector<int>& speciesPdg);

  void finalize();
  void validate(bool requireFinalized = true) const;

  int activityClass(double nch) const;
  int selectedMultiplicityClass(double nSelected) const;

  int numberOfSpecies() const;

  int numberOfCompositionPairs() const;
  std::uint64_t numberOfCompositionTemplates() const;
  int compositionPairIndex(int nch, int nSelected) const;

  TuneSpecies* speciesAt(int index);
  const TuneSpecies* speciesAt(int index) const;

  TuneSpecies* findSpecies(int pdg);
  const TuneSpecies* findSpecies(int pdg) const;

  void detachFromDirectories();

  void save(const std::string& fileName, const std::string& objectName = kObjectName) const;

  static std::unique_ptr<Tune> load(const std::string& fileName, const std::string& objectName = kObjectName);

  ClassDefOverride(Tune, 8);
};

} // namespace Ditto

#endif // Ditto_TUNE_H
