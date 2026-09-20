///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTune.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/14
/// \brief  Implementation of the Ditto tune data model.
///

#include "DittoTune.h"

#include <TFile.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <new>
#include <stdexcept>

ClassImp(Ditto::TuneSpecies);
ClassImp(Ditto::Tune);

namespace Ditto
{

namespace
{

void normalize1D(TH1D& histogram)
{
  const double integral = histogram.Integral(1, histogram.GetNbinsX());

  if (integral > 0.0) {
    histogram.Scale(1.0 / integral);
  }
}

void normalizeYSlices(TH2D& histogram)
{
  for (int ix = 1; ix <= histogram.GetNbinsX(); ++ix) {
    double sum = 0.0;

    for (int iy = 1; iy <= histogram.GetNbinsY(); ++iy) {
      sum += histogram.GetBinContent(ix, iy);
    }

    if (sum <= 0.0) {
      continue;
    }

    for (int iy = 1; iy <= histogram.GetNbinsY(); ++iy) {
      histogram.SetBinContent(ix, iy, histogram.GetBinContent(ix, iy) / sum);
      histogram.SetBinError(ix, iy, histogram.GetBinError(ix, iy) / sum);
    }
  }
}

void reconstructSparse(THnSparseD& target,
                       const char* name,
                       const char* title,
                       int nBinsX,
                       double xMin,
                       double xMax,
                       int nBinsY,
                       double yMin,
                       double yMax)
{
  const int bins[2] = {nBinsX, nBinsY};
  const double mins[2] = {xMin, yMin};
  const double maxs[2] = {xMax, yMax};

  target.~THnSparseD();
  new (&target) THnSparseD(name, title, 2, bins, mins, maxs);
}

void reconstructSparseLike(const THnSparseD& source,
                           THnSparseD& target)
{
  target.~THnSparseD();

  if (source.GetNdimensions() <= 0) {
    new (&target) THnSparseD();
    return;
  }

  const int nDimensions = source.GetNdimensions();
  std::vector<int> bins(static_cast<std::size_t>(nDimensions));
  std::vector<double> mins(static_cast<std::size_t>(nDimensions));
  std::vector<double> maxs(static_cast<std::size_t>(nDimensions));

  for (int i = 0; i < nDimensions; ++i) {
    const auto* axis = source.GetAxis(i);
    bins[static_cast<std::size_t>(i)] = axis->GetNbins();
    mins[static_cast<std::size_t>(i)] = axis->GetXmin();
    maxs[static_cast<std::size_t>(i)] = axis->GetXmax();
  }

  new (&target) THnSparseD(source.GetName(),
                           source.GetTitle(),
                           nDimensions,
                           bins.data(),
                           mins.data(),
                           maxs.data());

  for (int i = 0; i < nDimensions; ++i) {
    const auto* axis = source.GetAxis(i);
    const auto* variableBins = axis->GetXbins();

    if (variableBins && variableBins->GetSize() == axis->GetNbins() + 1) {
      target.SetBinEdges(i, variableBins->GetArray());
    }
  }

  target.Add(&source);
  target.SetEntries(source.GetEntries());
}

void copySparse(const THnSparseD& source,
                THnSparseD& target)
{
  target.Reset();
  target.Add(&source);
  target.SetEntries(source.GetEntries());
}

void copyNormalizedSparseYSlices(const THnSparseD& source,
                                 THnSparseD& target)
{
  if (source.GetNdimensions() != 2 ||
      target.GetNdimensions() != 2) {
    throw std::runtime_error("Ditto::TuneSpecies: multiplicity PMF must be two-dimensional");
  }

  target.Reset();

  const int nActivityBins = source.GetAxis(0)->GetNbins();
  const int nMultiplicityBins = source.GetAxis(1)->GetNbins();

  std::vector<double> sums(static_cast<std::size_t>(nActivityBins) + 2, 0.0);

  int coordinate[2] = {0, 0};

  const Long64_t nFilledBins = source.GetNbins();

  for (Long64_t i = 0; i < nFilledBins; ++i) {
    const double content = source.GetBinContent(i, coordinate);

    if (coordinate[0] < 1 ||
        coordinate[0] > nActivityBins ||
        coordinate[1] < 1 ||
        coordinate[1] > nMultiplicityBins) {
      continue;
    }

    sums[static_cast<std::size_t>(coordinate[0])] += content;
  }

  for (Long64_t i = 0; i < nFilledBins; ++i) {
    const double content = source.GetBinContent(i, coordinate);

    if (coordinate[0] < 1 ||
        coordinate[0] > nActivityBins ||
        coordinate[1] < 1 ||
        coordinate[1] > nMultiplicityBins) {
      continue;
    }

    const double sum = sums[static_cast<std::size_t>(coordinate[0])];

    if (sum <= 0.0 || content == 0.0)
      continue;

    const double point[2] = {
      source.GetAxis(0)->GetBinCenter(coordinate[0]),
      source.GetAxis(1)->GetBinCenter(coordinate[1])};

    target.Fill(point, content / sum);
  }

  target.SetEntries(source.GetEntries());
}

void finalizeMeanPt(const TH1D& ptSum,
                    const TH1D& particleCount,
                    TH1D& meanPt,
                    const std::string& name)
{
  meanPt = ptSum;
  meanPt.Reset("ICES");
  meanPt.SetName(name.c_str());

  for (int i = 1; i <= ptSum.GetNbinsX(); ++i) {
    const double count = particleCount.GetBinContent(i);

    if (count <= 0.0)
      continue;

    const double sumPt = ptSum.GetBinContent(i);
    const double mean = sumPt / count;

    meanPt.SetBinContent(i, mean);

    if (count <= 1.0)
      continue;

    // h*PtSumVsNch has Sumw2 enabled, therefore error^2 = sum(pT^2).
    const double sumPt2 =
      ptSum.GetBinError(i) *
      ptSum.GetBinError(i);

    const double variance =
      std::max(
        0.0,
        (sumPt2 - sumPt * sumPt / count) /
          (count - 1.0));

    meanPt.SetBinError(
      i,
      std::sqrt(variance / count));
  }
}

void normalizeSpeciesFractions(TObjArray& species,
                               bool centralCharged)
{
  if (species.GetEntriesFast() <= 0)
    return;

  auto* first = static_cast<TuneSpecies*>(species.UncheckedAt(0));

  if (!first)
    return;

  auto& reference = centralCharged
                      ? first->mPCentralChargedSpeciesGivenActivitySelected
                      : first->mPOtherSpeciesGivenActivitySelected;

  for (int ix = 1; ix <= reference.GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= reference.GetNbinsY(); ++iy) {
      double sum = 0.0;

      for (int iSpecies = 0; iSpecies < species.GetEntriesFast(); ++iSpecies) {
        auto* entry = static_cast<TuneSpecies*>(species.UncheckedAt(iSpecies));

        if (!entry)
          continue;

        const auto& histogram = centralCharged
                                  ? entry->mPCentralChargedSpeciesGivenActivitySelected
                                  : entry->mPOtherSpeciesGivenActivitySelected;

        sum += histogram.GetBinContent(ix, iy);
      }

      if (sum <= 0.0)
        continue;

      for (int iSpecies = 0; iSpecies < species.GetEntriesFast(); ++iSpecies) {
        auto* entry = static_cast<TuneSpecies*>(species.UncheckedAt(iSpecies));

        if (!entry)
          continue;

        auto& histogram = centralCharged
                            ? entry->mPCentralChargedSpeciesGivenActivitySelected
                            : entry->mPOtherSpeciesGivenActivitySelected;

        histogram.SetBinContent(ix,
                                iy,
                                histogram.GetBinContent(ix, iy) / sum);
        histogram.SetBinError(ix, iy, 0.0);
      }
    }
  }
}

} // namespace

TuneSpecies::TuneSpecies()
{
  detachFromDirectories();
}

TuneSpecies::TuneSpecies(int pdgIn,
                         int nActivityBins,
                         const double* activityEdges,
                         int nSelectedBins,
                         const double* selectedMultiplicityEdges,
                         int maxSpeciesMultiplicity,
                         int maxNch,
                         int nPtBins,
                         double ptMax,
                         int nEtaBins,
                         double etaMax)
  : mPdg(pdgIn),
    mHCountVsActivity(("hCountVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
                      ";N_{ch};N_{species}",
                      nActivityBins,
                      activityEdges,
                      maxSpeciesMultiplicity + 1,
                      -0.5,
                      maxSpeciesMultiplicity + 0.5),
    mHPtVsActivity(("hPtVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
                   ";N_{ch};p_{T} (GeV/c)",
                   nActivityBins,
                   activityEdges,
                   nPtBins,
                   0.0,
                   ptMax),
    mHEtaVsActivity(("hEtaVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
                    ";N_{ch};#eta",
                    nActivityBins,
                    activityEdges,
                    nEtaBins,
                    -etaMax,
                    etaMax),
    mHCentralChargedPtSumVsNch(("hCentralChargedPtSumVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                               ";N_{ch};#Sigma p_{T} (GeV/c)",
                               maxNch + 1,
                               -0.5,
                               maxNch + 0.5),
    mHCentralChargedPtCountVsNch(("hCentralChargedPtCountVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                                 ";N_{ch};particles",
                                 maxNch + 1,
                                 -0.5,
                                 maxNch + 0.5),
    mHOtherPtSumVsNch(("hOtherPtSumVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                      ";N_{ch};#Sigma p_{T} (GeV/c)",
                      maxNch + 1,
                      -0.5,
                      maxNch + 0.5),
    mHOtherPtCountVsNch(("hOtherPtCountVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                        ";N_{ch};particles",
                        maxNch + 1,
                        -0.5,
                        maxNch + 0.5),
    mHCentralChargedCountVsActivitySelected(("hCentralChargedCountVsActivitySelected_pdg_" + std::to_string(pdgIn)).c_str(),
                                            ";N_{ch};N_{selected}",
                                            nActivityBins,
                                            activityEdges,
                                            nSelectedBins,
                                            selectedMultiplicityEdges),
    mHOtherCountVsActivitySelected(("hOtherCountVsActivitySelected_pdg_" + std::to_string(pdgIn)).c_str(),
                                   ";N_{ch};N_{selected}",
                                   nActivityBins,
                                   activityEdges,
                                   nSelectedBins,
                                   selectedMultiplicityEdges),
    mHCentralChargedMultiplicityVsActivity(
      ("hCentralChargedMultiplicityVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{central}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data()),
    mHOtherMultiplicityVsActivity(
      ("hOtherMultiplicityVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{other}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data()),
    mPCentralChargedMultiplicityGivenActivity(
      ("pCentralChargedMultiplicityGivenActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{central}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data()),
    mPOtherMultiplicityGivenActivity(
      ("pOtherMultiplicityGivenActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{other}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data())
{
  mHCentralChargedMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  mHOtherMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  mPCentralChargedMultiplicityGivenActivity.SetBinEdges(0, activityEdges);
  mPOtherMultiplicityGivenActivity.SetBinEdges(0, activityEdges);

  mHCountVsActivity.Sumw2();
  mHPtVsActivity.Sumw2();
  mHEtaVsActivity.Sumw2();

  // For the pT-sum histograms Sumw2 stores sum(pT^2), which is used during
  // finalization to compute the statistical uncertainty on <pT>.
  mHCentralChargedPtSumVsNch.Sumw2();
  mHOtherPtSumVsNch.Sumw2();

  mHCentralChargedMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  mHOtherMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  mPCentralChargedMultiplicityGivenActivity.SetBinEdges(0, activityEdges);
  mPOtherMultiplicityGivenActivity.SetBinEdges(0, activityEdges);

  mHCentralChargedCountVsActivitySelected.Sumw2();
  mHOtherCountVsActivitySelected.Sumw2();

  // The sparse multiplicity tables are filled with unit event weights and are
  // used only as PMFs, so Sumw2 is intentionally not allocated for them.
  detachFromDirectories();
}

void TuneSpecies::finalize()
{
  mPCountGivenActivity = mHCountVsActivity;
  mPPtGivenActivity = mHPtVsActivity;
  mPEtaGivenActivity = mHEtaVsActivity;

  finalizeMeanPt(
    mHCentralChargedPtSumVsNch,
    mHCentralChargedPtCountVsNch,
    mHCentralChargedMeanPtVsNch,
    "hCentralChargedMeanPtVsNch_pdg_" + std::to_string(mPdg));

  finalizeMeanPt(
    mHOtherPtSumVsNch,
    mHOtherPtCountVsNch,
    mHOtherMeanPtVsNch,
    "hOtherMeanPtVsNch_pdg_" + std::to_string(mPdg));

  mPCentralChargedSpeciesGivenActivitySelected = mHCentralChargedCountVsActivitySelected;
  mPOtherSpeciesGivenActivitySelected = mHOtherCountVsActivitySelected;

  mPCountGivenActivity.SetName(("pCountGivenActivity_pdg_" + std::to_string(mPdg)).c_str());
  mPPtGivenActivity.SetName(("pPtGivenActivity_pdg_" + std::to_string(mPdg)).c_str());
  mPEtaGivenActivity.SetName(("pEtaGivenActivity_pdg_" + std::to_string(mPdg)).c_str());

  mPCentralChargedSpeciesGivenActivitySelected.SetName(("pCentralChargedSpeciesGivenActivitySelected_pdg_" + std::to_string(mPdg)).c_str());
  mPOtherSpeciesGivenActivitySelected.SetName(("pOtherSpeciesGivenActivitySelected_pdg_" + std::to_string(mPdg)).c_str());

  normalizeYSlices(mPCountGivenActivity);
  normalizeYSlices(mPPtGivenActivity);
  normalizeYSlices(mPEtaGivenActivity);

  copyNormalizedSparseYSlices(mHCentralChargedMultiplicityVsActivity,
                              mPCentralChargedMultiplicityGivenActivity);
  copyNormalizedSparseYSlices(mHOtherMultiplicityVsActivity,
                              mPOtherMultiplicityGivenActivity);

  detachFromDirectories();
}

void TuneSpecies::detachFromDirectories()
{
  mHCountVsActivity.SetDirectory(nullptr);
  mHPtVsActivity.SetDirectory(nullptr);
  mHEtaVsActivity.SetDirectory(nullptr);

  mHCentralChargedPtSumVsNch.SetDirectory(nullptr);
  mHCentralChargedPtCountVsNch.SetDirectory(nullptr);
  mHOtherPtSumVsNch.SetDirectory(nullptr);
  mHOtherPtCountVsNch.SetDirectory(nullptr);

  mHCentralChargedMeanPtVsNch.SetDirectory(nullptr);
  mHOtherMeanPtVsNch.SetDirectory(nullptr);

  mPCountGivenActivity.SetDirectory(nullptr);
  mPPtGivenActivity.SetDirectory(nullptr);
  mPEtaGivenActivity.SetDirectory(nullptr);

  mHCentralChargedCountVsActivitySelected.SetDirectory(nullptr);
  mHOtherCountVsActivitySelected.SetDirectory(nullptr);

  mPCentralChargedSpeciesGivenActivitySelected.SetDirectory(nullptr);
  mPOtherSpeciesGivenActivitySelected.SetDirectory(nullptr);
}

namespace
{

std::vector<double> histogramAxisEdges(const TAxis& axis)
{
  std::vector<double> edges(static_cast<std::size_t>(axis.GetNbins()) + 1);

  for (int i = 1; i <= axis.GetNbins(); ++i) {
    edges[static_cast<std::size_t>(i - 1)] = axis.GetBinLowEdge(i);
  }

  edges.back() = axis.GetBinUpEdge(axis.GetNbins());

  return edges;
}

TuneSpecies* cloneTuneSpecies(const TuneSpecies& source)
{
  const auto activityEdges =
    histogramAxisEdges(*source.mHCountVsActivity.GetXaxis());

  const auto selectedMultiplicityEdges =
    histogramAxisEdges(*source.mHCentralChargedCountVsActivitySelected.GetYaxis());

  auto* result = new TuneSpecies(
    source.mPdg,
    source.mHCountVsActivity.GetNbinsX(),
    activityEdges.data(),
    source.mHCentralChargedCountVsActivitySelected.GetNbinsY(),
    selectedMultiplicityEdges.data(),
    source.mHCountVsActivity.GetNbinsY() - 1,
    source.mHCentralChargedPtCountVsNch.GetNbinsX() - 1,
    source.mHPtVsActivity.GetNbinsY(),
    source.mHPtVsActivity.GetYaxis()->GetXmax(),
    source.mHEtaVsActivity.GetNbinsY(),
    source.mHEtaVsActivity.GetYaxis()->GetXmax());

  result->mParticleName = source.mParticleName;
  result->mMass = source.mMass;
  result->mChargeType = source.mChargeType;

  result->mHCountVsActivity = source.mHCountVsActivity;
  result->mHPtVsActivity = source.mHPtVsActivity;
  result->mHEtaVsActivity = source.mHEtaVsActivity;

  result->mPCountGivenActivity = source.mPCountGivenActivity;
  result->mPPtGivenActivity = source.mPPtGivenActivity;
  result->mPEtaGivenActivity = source.mPEtaGivenActivity;

  result->mHCentralChargedPtSumVsNch = source.mHCentralChargedPtSumVsNch;
  result->mHCentralChargedPtCountVsNch = source.mHCentralChargedPtCountVsNch;
  result->mHOtherPtSumVsNch = source.mHOtherPtSumVsNch;
  result->mHOtherPtCountVsNch = source.mHOtherPtCountVsNch;

  result->mHCentralChargedMeanPtVsNch = source.mHCentralChargedMeanPtVsNch;
  result->mHOtherMeanPtVsNch = source.mHOtherMeanPtVsNch;

  result->mHCentralChargedCountVsActivitySelected = source.mHCentralChargedCountVsActivitySelected;
  result->mHOtherCountVsActivitySelected = source.mHOtherCountVsActivitySelected;

  result->mPCentralChargedSpeciesGivenActivitySelected = source.mPCentralChargedSpeciesGivenActivitySelected;
  result->mPOtherSpeciesGivenActivitySelected = source.mPOtherSpeciesGivenActivitySelected;

  copySparse(source.mHCentralChargedMultiplicityVsActivity,
             result->mHCentralChargedMultiplicityVsActivity);
  copySparse(source.mHOtherMultiplicityVsActivity,
             result->mHOtherMultiplicityVsActivity);

  copySparse(source.mPCentralChargedMultiplicityGivenActivity,
             result->mPCentralChargedMultiplicityGivenActivity);
  copySparse(source.mPOtherMultiplicityGivenActivity,
             result->mPOtherMultiplicityGivenActivity);

  result->detachFromDirectories();

  return result;
}

} // namespace

Tune::Tune()
{
  mSpecies.SetOwner(kTRUE);
  detachFromDirectories();
}

Tune::Tune(const Tune& other)
  : TObject(other),
    mFormatVersion(other.mFormatVersion),
    mFinalized(other.mFinalized),
    mTeacher(other.mTeacher),
    mPythiaCard(other.mPythiaCard),
    mPythiaCardContent(other.mPythiaCardContent),
    mBeamIdA(other.mBeamIdA),
    mBeamIdB(other.mBeamIdB),
    mBeamFrameType(other.mBeamFrameType),
    mSqrtSNN(other.mSqrtSNN),
    mAzimuthModel(other.mAzimuthModel),
    mFinalStatus(other.mFinalStatus),
    mNEvents(other.mNEvents),
    mNAttempts(other.mNAttempts),
    mGenerationTimeSeconds(other.mGenerationTimeSeconds),
    mActivityOverflowEvents(other.mActivityOverflowEvents),
    mSelectedMultiplicityOverflowEvents(other.mSelectedMultiplicityOverflowEvents),
    mPtOverflowParticles(other.mPtOverflowParticles),
    mSpeciesMultiplicityOverflowEvents(other.mSpeciesMultiplicityOverflowEvents),
    mCentralChargedCoverageMismatchEvents(other.mCentralChargedCoverageMismatchEvents),
    mCentralChargedCoverageMissingParticles(other.mCentralChargedCoverageMissingParticles),
    mActivityEtaMax(other.mActivityEtaMax),
    mParticleEtaMax(other.mParticleEtaMax),
    mPtMax(other.mPtMax),
    mActivityEdges(other.mActivityEdges),
    mSelectedMultiplicityEdges(other.mSelectedMultiplicityEdges),
    mHNch(other.mHNch),
    mHNSelected(other.mHNSelected),
    mHNchSelected(other.mHNchSelected),
    mHEventsVsActivitySelected(other.mHEventsVsActivitySelected),
    mPNch(other.mPNch),
    mPNSelected(other.mPNSelected),
    mPNchSelected(other.mPNchSelected),
    mCompositionTemplateCapPerPair(other.mCompositionTemplateCapPerPair),
    mCompositionPairNch(other.mCompositionPairNch),
    mCompositionPairNSelected(other.mCompositionPairNSelected),
    mCompositionPairOffsets(other.mCompositionPairOffsets),
    mCompositionPairEventsSeen(other.mCompositionPairEventsSeen),
    mCompositionCentralCounts(other.mCompositionCentralCounts),
    mCompositionOtherCounts(other.mCompositionOtherCounts)
{
  reconstructSparseLike(other.mHNSelectedVsNch, mHNSelectedVsNch);
  reconstructSparseLike(other.mHNchSelectedVsNch, mHNchSelectedVsNch);
  reconstructSparseLike(other.mPNSelectedGivenNch, mPNSelectedGivenNch);
  reconstructSparseLike(other.mPNchSelectedGivenNch, mPNchSelectedGivenNch);

  mSpecies.SetOwner(kTRUE);

  for (int i = 0; i < other.numberOfSpecies(); ++i) {
    mSpecies.Add(cloneTuneSpecies(*other.speciesAt(i)));
  }

  detachFromDirectories();
}

Tune& Tune::operator=(const Tune& other)
{
  if (this == &other) {
    return *this;
  }

  TObject::operator=(other);

  mFormatVersion = other.mFormatVersion;
  mFinalized = other.mFinalized;

  mTeacher = other.mTeacher;
  mPythiaCard = other.mPythiaCard;
  mPythiaCardContent = other.mPythiaCardContent;

  mBeamIdA = other.mBeamIdA;
  mBeamIdB = other.mBeamIdB;
  mBeamFrameType = other.mBeamFrameType;
  mSqrtSNN = other.mSqrtSNN;
  mAzimuthModel = other.mAzimuthModel;
  mFinalStatus = other.mFinalStatus;

  mNEvents = other.mNEvents;
  mNAttempts = other.mNAttempts;
  mGenerationTimeSeconds = other.mGenerationTimeSeconds;

  mActivityOverflowEvents = other.mActivityOverflowEvents;
  mSelectedMultiplicityOverflowEvents = other.mSelectedMultiplicityOverflowEvents;
  mPtOverflowParticles = other.mPtOverflowParticles;
  mSpeciesMultiplicityOverflowEvents = other.mSpeciesMultiplicityOverflowEvents;

  mCentralChargedCoverageMismatchEvents = other.mCentralChargedCoverageMismatchEvents;
  mCentralChargedCoverageMissingParticles = other.mCentralChargedCoverageMissingParticles;

  mActivityEtaMax = other.mActivityEtaMax;
  mParticleEtaMax = other.mParticleEtaMax;
  mPtMax = other.mPtMax;

  mActivityEdges = other.mActivityEdges;
  mSelectedMultiplicityEdges = other.mSelectedMultiplicityEdges;

  mHNch = other.mHNch;
  mHNSelected = other.mHNSelected;
  mHNchSelected = other.mHNchSelected;

  reconstructSparseLike(other.mHNSelectedVsNch, mHNSelectedVsNch);
  reconstructSparseLike(other.mHNchSelectedVsNch, mHNchSelectedVsNch);
  mHEventsVsActivitySelected = other.mHEventsVsActivitySelected;

  mPNch = other.mPNch;
  mPNSelected = other.mPNSelected;
  mPNchSelected = other.mPNchSelected;

  reconstructSparseLike(other.mPNSelectedGivenNch, mPNSelectedGivenNch);
  reconstructSparseLike(other.mPNchSelectedGivenNch, mPNchSelectedGivenNch);

  mCompositionTemplateCapPerPair = other.mCompositionTemplateCapPerPair;
  mCompositionPairNch = other.mCompositionPairNch;
  mCompositionPairNSelected = other.mCompositionPairNSelected;
  mCompositionPairOffsets = other.mCompositionPairOffsets;
  mCompositionPairEventsSeen = other.mCompositionPairEventsSeen;
  mCompositionCentralCounts = other.mCompositionCentralCounts;
  mCompositionOtherCounts = other.mCompositionOtherCounts;

  mSpecies.Delete();
  mSpecies.SetOwner(kTRUE);

  for (int i = 0; i < other.numberOfSpecies(); ++i) {
    mSpecies.Add(cloneTuneSpecies(*other.speciesAt(i)));
  }

  detachFromDirectories();

  return *this;
}

Tune::~Tune()
{
  mSpecies.Delete();
}

void Tune::initialize(double activityEtaMaxIn,
                      double particleEtaMaxIn,
                      const std::vector<double>& activityEdgesIn,
                      const std::vector<double>& selectedMultiplicityEdgesIn,
                      int maxNch,
                      int maxSelectedMultiplicity,
                      int maxSpeciesMultiplicity,
                      int nPtBins,
                      double ptMaxIn,
                      int nEtaBins,
                      const std::vector<int>& speciesPdg)
{
  mFormatVersion = kFormatVersion;
  mFinalized = false;

  mTeacher = "PYTHIA8";
  mPythiaCard.clear();
  mPythiaCardContent.clear();

  mBeamIdA = 0;
  mBeamIdB = 0;
  mBeamFrameType = 0;
  mSqrtSNN = 0.0;
  mAzimuthModel = "uniform";
  mFinalStatus = 1;

  mActivityEtaMax = activityEtaMaxIn;
  mParticleEtaMax = particleEtaMaxIn;
  mPtMax = ptMaxIn;
  mActivityEdges = activityEdgesIn;
  mSelectedMultiplicityEdges = selectedMultiplicityEdgesIn;

  mHNch = TH1D("hNch", ";N_{ch};events", maxNch + 1, -0.5, maxNch + 0.5);

  mHNSelected = TH1D("hNSelected", ";N_{selected};events", maxSelectedMultiplicity + 1, -0.5, maxSelectedMultiplicity + 0.5);

  mHNchSelected = TH1D("hNchSelected", ";N_{ch}^{selected};events", maxNch + 1, -0.5, maxNch + 0.5);

  reconstructSparse(mHNSelectedVsNch,
                    "hNSelectedVsNch",
                    ";N_{ch};N_{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxSelectedMultiplicity + 1,
                    -0.5,
                    maxSelectedMultiplicity + 0.5);

  reconstructSparse(mHNchSelectedVsNch,
                    "hNchSelectedVsNch",
                    ";N_{ch};N_{ch}^{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5);

  const int nActivityBins = static_cast<int>(mActivityEdges.size()) - 1;
  const int nSelectedBins = static_cast<int>(mSelectedMultiplicityEdges.size()) - 1;

  mHEventsVsActivitySelected = TH2D("hEventsVsActivitySelected", ";N_{ch};N_{selected}", nActivityBins, mActivityEdges.data(), nSelectedBins, mSelectedMultiplicityEdges.data());

  mHNch.Sumw2();
  mHNSelected.Sumw2();
  mHNchSelected.Sumw2();

  mPNch = TH1D();
  mPNSelected = TH1D();
  mPNchSelected = TH1D();

  reconstructSparse(mPNSelectedGivenNch,
                    "pNSelectedGivenNch",
                    ";N_{ch};N_{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxSelectedMultiplicity + 1,
                    -0.5,
                    maxSelectedMultiplicity + 0.5);

  reconstructSparse(mPNchSelectedGivenNch,
                    "pNchSelectedGivenNch",
                    ";N_{ch};N_{ch}^{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5);

  mCompositionTemplateCapPerPair = 0;
  mCompositionPairNch.clear();
  mCompositionPairNSelected.clear();
  mCompositionPairOffsets.clear();
  mCompositionPairEventsSeen.clear();
  mCompositionCentralCounts.clear();
  mCompositionOtherCounts.clear();

  mSpecies.Delete();
  mSpecies.SetOwner(kTRUE);

  for (const int pdg : speciesPdg) {
    mSpecies.Add(new TuneSpecies(pdg,
                                 nActivityBins,
                                 mActivityEdges.data(),
                                 nSelectedBins,
                                 mSelectedMultiplicityEdges.data(),
                                 maxSpeciesMultiplicity,
                                 maxNch,
                                 nPtBins,
                                 mPtMax,
                                 nEtaBins,
                                 mParticleEtaMax));
  }

  detachFromDirectories();
}

void Tune::finalize()
{
  mPNch = mHNch;
  mPNSelected = mHNSelected;
  mPNchSelected = mHNchSelected;

  mPNch.SetName("pNch");
  mPNSelected.SetName("pNSelected");
  mPNchSelected.SetName("pNchSelected");

  normalize1D(mPNch);
  normalize1D(mPNSelected);
  normalize1D(mPNchSelected);

  copyNormalizedSparseYSlices(mHNSelectedVsNch, mPNSelectedGivenNch);
  copyNormalizedSparseYSlices(mHNchSelectedVsNch, mPNchSelectedGivenNch);

  for (int i = 0; i < numberOfSpecies(); ++i) {
    speciesAt(i)->finalize();
  }

  normalizeSpeciesFractions(mSpecies, true);
  normalizeSpeciesFractions(mSpecies, false);

  mFinalized = true;

  detachFromDirectories();
}

void Tune::validate(bool requireFinalized) const
{
  if (mFormatVersion != kFormatVersion) {
    throw std::runtime_error("Ditto::Tune: unsupported tune format version");
  }

  if (mTeacher.empty() ||
      mPythiaCardContent.empty()) {
    throw std::runtime_error("Ditto::Tune: missing teacher generator-card metadata");
  }

  if (mBeamIdA == 0 || mBeamIdB == 0) {
    throw std::runtime_error("Ditto::Tune: invalid beam IDs");
  }

  if (mBeamFrameType != 1) {
    throw std::runtime_error("Ditto::Tune: only Beams:frameType = 1 is supported");
  }

  if (!std::isfinite(mSqrtSNN) || mSqrtSNN <= 0.0) {
    throw std::runtime_error("Ditto::Tune: invalid sqrtSNN");
  }

  if (mAzimuthModel != "uniform") {
    throw std::runtime_error("Ditto::Tune: unsupported azimuth model");
  }

  if (mFinalStatus <= 0) {
    throw std::runtime_error("Ditto::Tune: finalStatus must be positive");
  }

  if (mActivityEtaMax <= 0.0) {
    throw std::runtime_error("Ditto::Tune: invalid activityEtaMax");
  }

  if (mParticleEtaMax <= 0.0) {
    throw std::runtime_error("Ditto::Tune: invalid particleEtaMax");
  }

  if (mActivityEdges.size() < 2) {
    throw std::runtime_error("Ditto::Tune: invalid activity edges");
  }

  if (mSelectedMultiplicityEdges.size() < 2) {
    throw std::runtime_error("Ditto::Tune: invalid selected-multiplicity edges");
  }

  for (std::size_t i = 1; i < mActivityEdges.size(); ++i) {
    if (mActivityEdges[i] <= mActivityEdges[i - 1]) {
      throw std::runtime_error("Ditto::Tune: activity edges are not strictly increasing");
    }
  }

  for (std::size_t i = 1; i < mSelectedMultiplicityEdges.size(); ++i) {
    if (mSelectedMultiplicityEdges[i] <= mSelectedMultiplicityEdges[i - 1]) {
      throw std::runtime_error("Ditto::Tune: selected-multiplicity edges are not strictly increasing");
    }
  }

  if (numberOfSpecies() <= 0) {
    throw std::runtime_error("Ditto::Tune: no species stored");
  }

  if (requireFinalized && !mFinalized) {
    throw std::runtime_error("Ditto::Tune: tune has not been finalized");
  }

  const int nActivityBins = static_cast<int>(mActivityEdges.size()) - 1;
  const int nSelectedBins = static_cast<int>(mSelectedMultiplicityEdges.size()) - 1;

  if (mHEventsVsActivitySelected.GetNbinsX() != nActivityBins ||
      mHEventsVsActivitySelected.GetNbinsY() != nSelectedBins) {
    throw std::runtime_error("Ditto::Tune: inconsistent event-composition binning");
  }

  if (requireFinalized) {
    if (mPNSelectedGivenNch.GetNdimensions() != 2 ||
        mPNchSelectedGivenNch.GetNdimensions() != 2 ||
        mPNSelectedGivenNch.GetAxis(0)->GetNbins() != mHNch.GetNbinsX() ||
        mPNSelectedGivenNch.GetAxis(1)->GetNbins() != mHNSelected.GetNbinsX() ||
        mPNchSelectedGivenNch.GetAxis(0)->GetNbins() != mHNch.GetNbinsX() ||
        mPNchSelectedGivenNch.GetAxis(1)->GetNbins() != mHNchSelected.GetNbinsX()) {
      throw std::runtime_error("Ditto::Tune: inconsistent sparse exact-multiplicity binning");
    }

    if (mCompositionPairNch.size() != mCompositionPairNSelected.size() ||
        mCompositionPairNch.size() != mCompositionPairEventsSeen.size()) {
      throw std::runtime_error("Ditto::Tune: inconsistent composition-pair metadata");
    }

    if (mCompositionPairOffsets.size() != mCompositionPairNch.size() + 1 ||
        mCompositionPairOffsets.empty() ||
        mCompositionPairOffsets.front() != 0) {
      throw std::runtime_error("Ditto::Tune: invalid composition-template offsets");
    }

    for (std::size_t i = 1; i < mCompositionPairOffsets.size(); ++i) {
      if (mCompositionPairOffsets[i] < mCompositionPairOffsets[i - 1]) {
        throw std::runtime_error("Ditto::Tune: composition-template offsets are not monotonic");
      }
    }

    for (std::size_t i = 1; i < mCompositionPairNch.size(); ++i) {
      if (mCompositionPairNch[i] < mCompositionPairNch[i - 1] ||
          (mCompositionPairNch[i] == mCompositionPairNch[i - 1] &&
           mCompositionPairNSelected[i] <= mCompositionPairNSelected[i - 1])) {
        throw std::runtime_error("Ditto::Tune: composition pairs are not strictly ordered");
      }
    }

    if (mCompositionTemplateCapPerPair == 0) {
      throw std::runtime_error("Ditto::Tune: invalid composition-template cap");
    }

    const std::uint64_t nTemplates =
      mCompositionPairOffsets.back();

    const std::uint64_t expectedCounts =
      nTemplates *
      static_cast<std::uint64_t>(numberOfSpecies());

    if (mCompositionCentralCounts.size() != expectedCounts ||
        mCompositionOtherCounts.size() != expectedCounts) {
      throw std::runtime_error("Ditto::Tune: inconsistent flattened composition-template size");
    }

    const std::size_t nSpecies =
      static_cast<std::size_t>(numberOfSpecies());

    for (std::size_t iPair = 0;
         iPair < mCompositionPairNch.size();
         ++iPair) {
      const std::uint64_t first =
        mCompositionPairOffsets[iPair];

      const std::uint64_t last =
        mCompositionPairOffsets[iPair + 1];

      const std::uint64_t stored =
        last - first;

      if (stored == 0 ||
          stored > mCompositionTemplateCapPerPair ||
          stored > mCompositionPairEventsSeen[iPair]) {
        throw std::runtime_error("Ditto::Tune: invalid composition-template reservoir metadata");
      }

      for (std::uint64_t iTemplate = first;
           iTemplate < last;
           ++iTemplate) {
        const std::uint64_t base =
          iTemplate *
          static_cast<std::uint64_t>(nSpecies);

        int centralTotal = 0;
        int selectedTotal = 0;

        for (std::size_t iSpecies = 0;
             iSpecies < nSpecies;
             ++iSpecies) {
          const int central =
            static_cast<int>(
              mCompositionCentralCounts[static_cast<std::size_t>(base) + iSpecies]);

          const int other =
            static_cast<int>(
              mCompositionOtherCounts[static_cast<std::size_t>(base) + iSpecies]);

          centralTotal += central;
          selectedTotal += central + other;
        }

        if (centralTotal > mCompositionPairNch[iPair] ||
            selectedTotal != mCompositionPairNSelected[iPair]) {
          throw std::runtime_error("Ditto::Tune: inconsistent empirical composition template");
        }
      }
    }
  }

  for (int i = 0; i < numberOfSpecies(); ++i) {
    const auto* entry = speciesAt(i);

    if (!entry) {
      throw std::runtime_error("Ditto::Tune: invalid species entry");
    }

    if (!std::isfinite(entry->mMass) ||
        entry->mMass < 0.0) {
      throw std::runtime_error("Ditto::Tune: invalid stored species mass");
    }

    if (entry->mHCountVsActivity.GetNbinsX() != nActivityBins ||
        entry->mHPtVsActivity.GetNbinsX() != nActivityBins ||
        entry->mHEtaVsActivity.GetNbinsX() != nActivityBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent activity binning");
    }

    const int nExactNchBins =
      mHNch.GetNbinsX();

    if (entry->mHCentralChargedPtSumVsNch.GetNbinsX() != nExactNchBins ||
        entry->mHCentralChargedPtCountVsNch.GetNbinsX() != nExactNchBins ||
        entry->mHOtherPtSumVsNch.GetNbinsX() != nExactNchBins ||
        entry->mHOtherPtCountVsNch.GetNbinsX() != nExactNchBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent exact-Nch pT-moment binning");
    }

    if (requireFinalized &&
        (entry->mHCentralChargedMeanPtVsNch.GetNbinsX() != nExactNchBins ||
         entry->mHOtherMeanPtVsNch.GetNbinsX() != nExactNchBins)) {
      throw std::runtime_error("Ditto::Tune: incomplete exact-Nch mean-pT tables");
    }

    if (entry->mHCentralChargedCountVsActivitySelected.GetNbinsX() != nActivityBins ||
        entry->mHCentralChargedCountVsActivitySelected.GetNbinsY() != nSelectedBins ||
        entry->mHOtherCountVsActivitySelected.GetNbinsX() != nActivityBins ||
        entry->mHOtherCountVsActivitySelected.GetNbinsY() != nSelectedBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent constrained-composition binning");
    }

    const int nSpeciesMultiplicityBins =
      entry->mHCountVsActivity.GetNbinsY();

    if (entry->mHCentralChargedMultiplicityVsActivity.GetNdimensions() != 2 ||
        entry->mHOtherMultiplicityVsActivity.GetNdimensions() != 2 ||
        entry->mHCentralChargedMultiplicityVsActivity.GetAxis(0)->GetNbins() != nActivityBins ||
        entry->mHOtherMultiplicityVsActivity.GetAxis(0)->GetNbins() != nActivityBins ||
        entry->mHCentralChargedMultiplicityVsActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins ||
        entry->mHOtherMultiplicityVsActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent sparse multiplicity-diagnostic binning");
    }

    if (requireFinalized &&
        (entry->mPCountGivenActivity.GetNbinsX() != nActivityBins ||
         entry->mPPtGivenActivity.GetNbinsX() != nActivityBins ||
         entry->mPEtaGivenActivity.GetNbinsX() != nActivityBins ||
         entry->mPCentralChargedSpeciesGivenActivitySelected.GetNbinsX() != nActivityBins ||
         entry->mPCentralChargedSpeciesGivenActivitySelected.GetNbinsY() != nSelectedBins ||
         entry->mPOtherSpeciesGivenActivitySelected.GetNbinsX() != nActivityBins ||
         entry->mPOtherSpeciesGivenActivitySelected.GetNbinsY() != nSelectedBins ||
         entry->mPCentralChargedMultiplicityGivenActivity.GetNdimensions() != 2 ||
         entry->mPOtherMultiplicityGivenActivity.GetNdimensions() != 2 ||
         entry->mPCentralChargedMultiplicityGivenActivity.GetAxis(0)->GetNbins() != nActivityBins ||
         entry->mPOtherMultiplicityGivenActivity.GetAxis(0)->GetNbins() != nActivityBins ||
         entry->mPCentralChargedMultiplicityGivenActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins ||
         entry->mPOtherMultiplicityGivenActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins)) {
      throw std::runtime_error("Ditto::Tune: incomplete probability tables");
    }
  }
}

int Tune::activityClass(double nch) const
{
  if (mActivityEdges.size() < 2 ||
      nch < mActivityEdges.front() ||
      nch >= mActivityEdges.back()) {
    return -1;
  }

  const auto upper = std::upper_bound(mActivityEdges.begin(), mActivityEdges.end(), nch);

  return static_cast<int>(std::distance(mActivityEdges.begin(), upper) - 1);
}

int Tune::selectedMultiplicityClass(double nSelected) const
{
  if (mSelectedMultiplicityEdges.size() < 2 ||
      nSelected < mSelectedMultiplicityEdges.front() ||
      nSelected >= mSelectedMultiplicityEdges.back()) {
    return -1;
  }

  const auto upper = std::upper_bound(mSelectedMultiplicityEdges.begin(),
                                      mSelectedMultiplicityEdges.end(),
                                      nSelected);

  return static_cast<int>(std::distance(mSelectedMultiplicityEdges.begin(), upper) - 1);
}

int Tune::numberOfSpecies() const
{
  return mSpecies.GetEntriesFast();
}

int Tune::numberOfCompositionPairs() const
{
  return static_cast<int>(mCompositionPairNch.size());
}

std::uint64_t Tune::numberOfCompositionTemplates() const
{
  return mCompositionPairOffsets.empty()
           ? 0
           : mCompositionPairOffsets.back();
}

int Tune::compositionPairIndex(int nch, int nSelected) const
{
  const auto lower =
    std::lower_bound(
      mCompositionPairNch.begin(),
      mCompositionPairNch.end(),
      nch);

  std::size_t index =
    static_cast<std::size_t>(
      std::distance(mCompositionPairNch.begin(), lower));

  while (index < mCompositionPairNch.size() &&
         mCompositionPairNch[index] == nch) {
    if (mCompositionPairNSelected[index] == nSelected) {
      return static_cast<int>(index);
    }

    if (mCompositionPairNSelected[index] > nSelected) {
      break;
    }

    ++index;
  }

  return -1;
}

TuneSpecies* Tune::speciesAt(int index)
{
  return static_cast<TuneSpecies*>(mSpecies.UncheckedAt(index));
}

const TuneSpecies* Tune::speciesAt(int index) const
{
  return static_cast<const TuneSpecies*>(mSpecies.UncheckedAt(index));
}

TuneSpecies* Tune::findSpecies(int pdg)
{
  for (int i = 0; i < numberOfSpecies(); ++i) {
    auto* entry = speciesAt(i);

    if (entry && entry->mPdg == pdg) {
      return entry;
    }
  }

  return nullptr;
}

const TuneSpecies* Tune::findSpecies(int pdg) const
{
  for (int i = 0; i < numberOfSpecies(); ++i) {
    const auto* entry = speciesAt(i);

    if (entry && entry->mPdg == pdg) {
      return entry;
    }
  }

  return nullptr;
}

void Tune::detachFromDirectories()
{
  mHNch.SetDirectory(nullptr);
  mHNSelected.SetDirectory(nullptr);
  mHNchSelected.SetDirectory(nullptr);

  mHEventsVsActivitySelected.SetDirectory(nullptr);

  mPNch.SetDirectory(nullptr);
  mPNSelected.SetDirectory(nullptr);
  mPNchSelected.SetDirectory(nullptr);

  mSpecies.SetOwner(kTRUE);

  for (int i = 0; i < numberOfSpecies(); ++i) {
    auto* entry = speciesAt(i);

    if (entry) {
      entry->detachFromDirectories();
    }
  }
}

void Tune::save(const std::string& fileName, const std::string& objectName) const
{
  validate(true);

  std::unique_ptr<TFile> output(TFile::Open(fileName.c_str(), "RECREATE"));

  if (!output || output->IsZombie()) {
    throw std::runtime_error("Ditto::Tune: could not create " + fileName);
  }

  output->WriteObject(this, objectName.c_str());
  output->Write();
  output->Close();
}

std::unique_ptr<Tune> Tune::load(const std::string& fileName,
                                 const std::string& objectName)
{
  std::unique_ptr<TFile> input(TFile::Open(fileName.c_str(), "READ"));

  if (!input || input->IsZombie()) {
    throw std::runtime_error("Ditto::Tune: could not open " + fileName);
  }

  Tune* stored = nullptr;
  input->GetObject(objectName.c_str(), stored);

  if (!stored) {
    throw std::runtime_error("Ditto::Tune: object '" + objectName + "' not found in " + fileName);
  }

  auto result = std::make_unique<Tune>(*stored);

  result->detachFromDirectories();
  result->validate(true);

  return result;
}

} // namespace Ditto
