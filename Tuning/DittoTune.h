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
  int mPdg = 0;
  std::string mParticleName;

  // Particle properties copied from the initialized PYTHIA teacher so the
  // Ditto generator does not need a runtime particle-data lookup.
  double mMass = 0.0;
  int mChargeType = 0; // PYTHIA convention: three times the electric charge.

  TH2D mHCountVsActivity;
  TH2D mHPtVsActivity;
  TH2D mHEtaVsActivity;

  TH2D mPCountGivenActivity;
  TH2D mPPtGivenActivity;
  TH2D mPEtaGivenActivity;

  // Exact-Nch, species-dependent pT moments. These retain the multiplicity
  // dependence washed out by the coarser activity classes used by
  // mPPtGivenActivity.
  //
  // "Central charged" and "other" use the same component definition as the
  // empirical composition templates. The pT-sum histograms have Sumw2 enabled:
  // their squared bin errors therefore store sum(pT^2), allowing the finalized
  // mean-pT histograms to carry the statistical error on the mean.
  TH1D mHCentralChargedPtSumVsNch;
  TH1D mHCentralChargedPtCountVsNch;
  TH1D mHOtherPtSumVsNch;
  TH1D mHOtherPtCountVsNch;

  TH1D mHCentralChargedMeanPtVsNch;
  TH1D mHOtherMeanPtVsNch;

  // Particle-count accumulators used to learn a constrained event
  // composition. "Central charged" means a charged selected particle inside
  // |eta| < mActivityEtaMax. "Other" contains every other selected particle.
  //
  // The axes are:
  //   x = Nch activity class
  //   y = Nselected class
  TH2D mHCentralChargedCountVsActivitySelected;
  TH2D mHOtherCountVsActivitySelected;

  // Species probabilities normalized across all species in each
  // (activity, Nselected) cell. Retained as useful composition diagnostics.
  TH2D mPCentralChargedSpeciesGivenActivitySelected;
  TH2D mPOtherSpeciesGivenActivitySelected;

  // Full per-species multiplicity distributions retained as diagnostics.
  // Generation uses the empirical correlated composition templates below.
  // These are filled once per event, including zero multiplicity. THnSparseD
  // keeps the large exact-multiplicity axis compact for heavy-ion tunes.
  //
  // Axes:
  //   0 = Nch activity class
  //   1 = exact species multiplicity
  THnSparseD mHCentralChargedMultiplicityVsActivity;
  THnSparseD mHOtherMultiplicityVsActivity;

  THnSparseD mPCentralChargedMultiplicityGivenActivity;
  THnSparseD mPOtherMultiplicityGivenActivity;

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

  int mFormatVersion = kFormatVersion;
  bool mFinalized = false;

  std::string mTeacher = "PYTHIA8";
  std::string mPythiaCard;
  std::string mPythiaCardContent;

  // Complete beam / generator metadata copied from the initialized teacher.
  // Ditto does not accept an independent collision system or energy: these
  // values are the generator configuration.
  int mBeamIdA = 0;
  int mBeamIdB = 0;
  int mBeamFrameType = 0;
  double mSqrtSNN = 0.0;

  // Current Ditto kinematic model samples the absolute azimuth uniformly.
  // Keeping this explicit in the card prevents an implicit analytic fallback.
  std::string mAzimuthModel = "uniform";

  // Final-state status written to TParticle / exported PYTHIA events.
  int mFinalStatus = 1;

  std::uint64_t mNEvents = 0;
  std::uint64_t mNAttempts = 0;
  double mGenerationTimeSeconds = 0.0;

  std::uint64_t mActivityOverflowEvents = 0;
  std::uint64_t mSelectedMultiplicityOverflowEvents = 0;
  std::uint64_t mPtOverflowParticles = 0;
  std::uint64_t mSpeciesMultiplicityOverflowEvents = 0;

  // Diagnostic for the exact Nch constraint. These count events/particles
  // for which final charged PYTHIA particles in the activity acceptance are
  // not represented by the selected species list.
  std::uint64_t mCentralChargedCoverageMismatchEvents = 0;
  std::uint64_t mCentralChargedCoverageMissingParticles = 0;

  double mActivityEtaMax = 0.0;
  double mParticleEtaMax = 0.0;
  double mPtMax = 0.0;

  std::vector<double> mActivityEdges;
  std::vector<double> mSelectedMultiplicityEdges;

  TH1D mHNch;
  TH1D mHNSelected;
  TH1D mHNchSelected;

  // Exact event-level correlations. PbPb occupies only a tiny fraction of the
  // full integer (Nch, Nselected) plane, so these are sparse. The raw
  // accumulators are tuner-only; the normalized conditional distributions
  // below are persisted and used by the generator.
  THnSparseD mHNSelectedVsNch;   //!
  THnSparseD mHNchSelectedVsNch; //!

  // Event occupancy for the coarser composition conditioning.
  TH2D mHEventsVsActivitySelected;

  TH1D mPNch;
  TH1D mPNSelected;
  TH1D mPNchSelected;

  THnSparseD mPNSelectedGivenNch;
  THnSparseD mPNchSelectedGivenNch;

  // Empirical correlated species-composition templates.
  //
  // Pairs are stored in lexicographic (Nch, Nselected) order. For pair i,
  // templates are in the half-open range
  //
  //   [mCompositionPairOffsets[i], mCompositionPairOffsets[i + 1])
  //
  // and each template contains numberOfSpecies() central and other counts in
  // the two flattened arrays below.
  std::uint64_t mCompositionTemplateCapPerPair = 0;

  std::vector<int> mCompositionPairNch;
  std::vector<int> mCompositionPairNSelected;
  std::vector<std::uint64_t> mCompositionPairOffsets;
  std::vector<std::uint64_t> mCompositionPairEventsSeen;

  std::vector<unsigned short> mCompositionCentralCounts;
  std::vector<unsigned short> mCompositionOtherCounts;

  // Owned TuneSpecies objects.
  TObjArray mSpecies;

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
