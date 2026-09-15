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
                      ? first->pCentralChargedSpeciesGivenActivitySelected
                      : first->pOtherSpeciesGivenActivitySelected;

  for (int ix = 1; ix <= reference.GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= reference.GetNbinsY(); ++iy) {
      double sum = 0.0;

      for (int iSpecies = 0; iSpecies < species.GetEntriesFast(); ++iSpecies) {
        auto* entry = static_cast<TuneSpecies*>(species.UncheckedAt(iSpecies));

        if (!entry)
          continue;

        const auto& histogram = centralCharged
                                  ? entry->pCentralChargedSpeciesGivenActivitySelected
                                  : entry->pOtherSpeciesGivenActivitySelected;

        sum += histogram.GetBinContent(ix, iy);
      }

      if (sum <= 0.0)
        continue;

      for (int iSpecies = 0; iSpecies < species.GetEntriesFast(); ++iSpecies) {
        auto* entry = static_cast<TuneSpecies*>(species.UncheckedAt(iSpecies));

        if (!entry)
          continue;

        auto& histogram = centralCharged
                            ? entry->pCentralChargedSpeciesGivenActivitySelected
                            : entry->pOtherSpeciesGivenActivitySelected;

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
  : pdg(pdgIn),
    hCountVsActivity(("hCountVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
                     ";N_{ch};N_{species}",
                     nActivityBins,
                     activityEdges,
                     maxSpeciesMultiplicity + 1,
                     -0.5,
                     maxSpeciesMultiplicity + 0.5),
    hPtVsActivity(("hPtVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
                  ";N_{ch};p_{T} (GeV/c)",
                  nActivityBins,
                  activityEdges,
                  nPtBins,
                  0.0,
                  ptMax),
    hEtaVsActivity(("hEtaVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
                   ";N_{ch};#eta",
                   nActivityBins,
                   activityEdges,
                   nEtaBins,
                   -etaMax,
                   etaMax),
    hCentralChargedPtSumVsNch(("hCentralChargedPtSumVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                              ";N_{ch};#Sigma p_{T} (GeV/c)",
                              maxNch + 1,
                              -0.5,
                              maxNch + 0.5),
    hCentralChargedPtCountVsNch(("hCentralChargedPtCountVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                                ";N_{ch};particles",
                                maxNch + 1,
                                -0.5,
                                maxNch + 0.5),
    hOtherPtSumVsNch(("hOtherPtSumVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                     ";N_{ch};#Sigma p_{T} (GeV/c)",
                     maxNch + 1,
                     -0.5,
                     maxNch + 0.5),
    hOtherPtCountVsNch(("hOtherPtCountVsNch_pdg_" + std::to_string(pdgIn)).c_str(),
                       ";N_{ch};particles",
                       maxNch + 1,
                       -0.5,
                       maxNch + 0.5),
    hCentralChargedCountVsActivitySelected(("hCentralChargedCountVsActivitySelected_pdg_" + std::to_string(pdgIn)).c_str(),
                                           ";N_{ch};N_{selected}",
                                           nActivityBins,
                                           activityEdges,
                                           nSelectedBins,
                                           selectedMultiplicityEdges),
    hOtherCountVsActivitySelected(("hOtherCountVsActivitySelected_pdg_" + std::to_string(pdgIn)).c_str(),
                                  ";N_{ch};N_{selected}",
                                  nActivityBins,
                                  activityEdges,
                                  nSelectedBins,
                                  selectedMultiplicityEdges),
    hCentralChargedMultiplicityVsActivity(
      ("hCentralChargedMultiplicityVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{central}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data()),
    hOtherMultiplicityVsActivity(
      ("hOtherMultiplicityVsActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{other}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data()),
    pCentralChargedMultiplicityGivenActivity(
      ("pCentralChargedMultiplicityGivenActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{central}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data()),
    pOtherMultiplicityGivenActivity(
      ("pOtherMultiplicityGivenActivity_pdg_" + std::to_string(pdgIn)).c_str(),
      ";N_{ch};N_{species}^{other}",
      2,
      std::array<int, 2>{nActivityBins, maxSpeciesMultiplicity + 1}.data(),
      std::array<double, 2>{activityEdges[0], -0.5}.data(),
      std::array<double, 2>{activityEdges[nActivityBins], maxSpeciesMultiplicity + 0.5}.data())
{
  hCentralChargedMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  hOtherMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  pCentralChargedMultiplicityGivenActivity.SetBinEdges(0, activityEdges);
  pOtherMultiplicityGivenActivity.SetBinEdges(0, activityEdges);

  hCountVsActivity.Sumw2();
  hPtVsActivity.Sumw2();
  hEtaVsActivity.Sumw2();

  // For the pT-sum histograms Sumw2 stores sum(pT^2), which is used during
  // finalization to compute the statistical uncertainty on <pT>.
  hCentralChargedPtSumVsNch.Sumw2();
  hOtherPtSumVsNch.Sumw2();

  hCentralChargedMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  hOtherMultiplicityVsActivity.SetBinEdges(0, activityEdges);
  pCentralChargedMultiplicityGivenActivity.SetBinEdges(0, activityEdges);
  pOtherMultiplicityGivenActivity.SetBinEdges(0, activityEdges);

  hCentralChargedCountVsActivitySelected.Sumw2();
  hOtherCountVsActivitySelected.Sumw2();

  // The sparse multiplicity tables are filled with unit event weights and are
  // used only as PMFs, so Sumw2 is intentionally not allocated for them.
  detachFromDirectories();
}

void TuneSpecies::finalize()
{
  pCountGivenActivity = hCountVsActivity;
  pPtGivenActivity = hPtVsActivity;
  pEtaGivenActivity = hEtaVsActivity;

  finalizeMeanPt(
    hCentralChargedPtSumVsNch,
    hCentralChargedPtCountVsNch,
    hCentralChargedMeanPtVsNch,
    "hCentralChargedMeanPtVsNch_pdg_" + std::to_string(pdg));

  finalizeMeanPt(
    hOtherPtSumVsNch,
    hOtherPtCountVsNch,
    hOtherMeanPtVsNch,
    "hOtherMeanPtVsNch_pdg_" + std::to_string(pdg));

  pCentralChargedSpeciesGivenActivitySelected = hCentralChargedCountVsActivitySelected;
  pOtherSpeciesGivenActivitySelected = hOtherCountVsActivitySelected;

  pCountGivenActivity.SetName(("pCountGivenActivity_pdg_" + std::to_string(pdg)).c_str());
  pPtGivenActivity.SetName(("pPtGivenActivity_pdg_" + std::to_string(pdg)).c_str());
  pEtaGivenActivity.SetName(("pEtaGivenActivity_pdg_" + std::to_string(pdg)).c_str());

  pCentralChargedSpeciesGivenActivitySelected.SetName(("pCentralChargedSpeciesGivenActivitySelected_pdg_" + std::to_string(pdg)).c_str());
  pOtherSpeciesGivenActivitySelected.SetName(("pOtherSpeciesGivenActivitySelected_pdg_" + std::to_string(pdg)).c_str());

  normalizeYSlices(pCountGivenActivity);
  normalizeYSlices(pPtGivenActivity);
  normalizeYSlices(pEtaGivenActivity);

  copyNormalizedSparseYSlices(hCentralChargedMultiplicityVsActivity,
                              pCentralChargedMultiplicityGivenActivity);
  copyNormalizedSparseYSlices(hOtherMultiplicityVsActivity,
                              pOtherMultiplicityGivenActivity);

  detachFromDirectories();
}

void TuneSpecies::detachFromDirectories()
{
  hCountVsActivity.SetDirectory(nullptr);
  hPtVsActivity.SetDirectory(nullptr);
  hEtaVsActivity.SetDirectory(nullptr);

  hCentralChargedPtSumVsNch.SetDirectory(nullptr);
  hCentralChargedPtCountVsNch.SetDirectory(nullptr);
  hOtherPtSumVsNch.SetDirectory(nullptr);
  hOtherPtCountVsNch.SetDirectory(nullptr);

  hCentralChargedMeanPtVsNch.SetDirectory(nullptr);
  hOtherMeanPtVsNch.SetDirectory(nullptr);

  pCountGivenActivity.SetDirectory(nullptr);
  pPtGivenActivity.SetDirectory(nullptr);
  pEtaGivenActivity.SetDirectory(nullptr);

  hCentralChargedCountVsActivitySelected.SetDirectory(nullptr);
  hOtherCountVsActivitySelected.SetDirectory(nullptr);

  pCentralChargedSpeciesGivenActivitySelected.SetDirectory(nullptr);
  pOtherSpeciesGivenActivitySelected.SetDirectory(nullptr);
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
    histogramAxisEdges(*source.hCountVsActivity.GetXaxis());

  const auto selectedMultiplicityEdges =
    histogramAxisEdges(*source.hCentralChargedCountVsActivitySelected.GetYaxis());

  auto* result = new TuneSpecies(
    source.pdg,
    source.hCountVsActivity.GetNbinsX(),
    activityEdges.data(),
    source.hCentralChargedCountVsActivitySelected.GetNbinsY(),
    selectedMultiplicityEdges.data(),
    source.hCountVsActivity.GetNbinsY() - 1,
    source.hCentralChargedPtCountVsNch.GetNbinsX() - 1,
    source.hPtVsActivity.GetNbinsY(),
    source.hPtVsActivity.GetYaxis()->GetXmax(),
    source.hEtaVsActivity.GetNbinsY(),
    source.hEtaVsActivity.GetYaxis()->GetXmax());

  result->particleName = source.particleName;
  result->mass = source.mass;
  result->chargeType = source.chargeType;

  result->hCountVsActivity = source.hCountVsActivity;
  result->hPtVsActivity = source.hPtVsActivity;
  result->hEtaVsActivity = source.hEtaVsActivity;

  result->pCountGivenActivity = source.pCountGivenActivity;
  result->pPtGivenActivity = source.pPtGivenActivity;
  result->pEtaGivenActivity = source.pEtaGivenActivity;

  result->hCentralChargedPtSumVsNch = source.hCentralChargedPtSumVsNch;
  result->hCentralChargedPtCountVsNch = source.hCentralChargedPtCountVsNch;
  result->hOtherPtSumVsNch = source.hOtherPtSumVsNch;
  result->hOtherPtCountVsNch = source.hOtherPtCountVsNch;

  result->hCentralChargedMeanPtVsNch = source.hCentralChargedMeanPtVsNch;
  result->hOtherMeanPtVsNch = source.hOtherMeanPtVsNch;

  result->hCentralChargedCountVsActivitySelected = source.hCentralChargedCountVsActivitySelected;
  result->hOtherCountVsActivitySelected = source.hOtherCountVsActivitySelected;

  result->pCentralChargedSpeciesGivenActivitySelected = source.pCentralChargedSpeciesGivenActivitySelected;
  result->pOtherSpeciesGivenActivitySelected = source.pOtherSpeciesGivenActivitySelected;

  copySparse(source.hCentralChargedMultiplicityVsActivity,
             result->hCentralChargedMultiplicityVsActivity);
  copySparse(source.hOtherMultiplicityVsActivity,
             result->hOtherMultiplicityVsActivity);

  copySparse(source.pCentralChargedMultiplicityGivenActivity,
             result->pCentralChargedMultiplicityGivenActivity);
  copySparse(source.pOtherMultiplicityGivenActivity,
             result->pOtherMultiplicityGivenActivity);

  result->detachFromDirectories();

  return result;
}

} // namespace

Tune::Tune()
{
  species.SetOwner(kTRUE);
  detachFromDirectories();
}

Tune::Tune(const Tune& other)
  : TObject(other),
    formatVersion(other.formatVersion),
    finalized(other.finalized),
    teacher(other.teacher),
    pythiaCard(other.pythiaCard),
    pythiaCardContent(other.pythiaCardContent),
    beamIdA(other.beamIdA),
    beamIdB(other.beamIdB),
    beamFrameType(other.beamFrameType),
    sqrtSNN(other.sqrtSNN),
    azimuthModel(other.azimuthModel),
    finalStatus(other.finalStatus),
    nEvents(other.nEvents),
    nAttempts(other.nAttempts),
    generationTimeSeconds(other.generationTimeSeconds),
    activityOverflowEvents(other.activityOverflowEvents),
    selectedMultiplicityOverflowEvents(other.selectedMultiplicityOverflowEvents),
    ptOverflowParticles(other.ptOverflowParticles),
    speciesMultiplicityOverflowEvents(other.speciesMultiplicityOverflowEvents),
    centralChargedCoverageMismatchEvents(other.centralChargedCoverageMismatchEvents),
    centralChargedCoverageMissingParticles(other.centralChargedCoverageMissingParticles),
    activityEtaMax(other.activityEtaMax),
    particleEtaMax(other.particleEtaMax),
    ptMax(other.ptMax),
    activityEdges(other.activityEdges),
    selectedMultiplicityEdges(other.selectedMultiplicityEdges),
    hNch(other.hNch),
    hNSelected(other.hNSelected),
    hNchSelected(other.hNchSelected),
    hEventsVsActivitySelected(other.hEventsVsActivitySelected),
    pNch(other.pNch),
    pNSelected(other.pNSelected),
    pNchSelected(other.pNchSelected),
    compositionTemplateCapPerPair(other.compositionTemplateCapPerPair),
    compositionPairNch(other.compositionPairNch),
    compositionPairNSelected(other.compositionPairNSelected),
    compositionPairOffsets(other.compositionPairOffsets),
    compositionPairEventsSeen(other.compositionPairEventsSeen),
    compositionCentralCounts(other.compositionCentralCounts),
    compositionOtherCounts(other.compositionOtherCounts)
{
  reconstructSparseLike(other.hNSelectedVsNch, hNSelectedVsNch);
  reconstructSparseLike(other.hNchSelectedVsNch, hNchSelectedVsNch);
  reconstructSparseLike(other.pNSelectedGivenNch, pNSelectedGivenNch);
  reconstructSparseLike(other.pNchSelectedGivenNch, pNchSelectedGivenNch);

  species.SetOwner(kTRUE);

  for (int i = 0; i < other.numberOfSpecies(); ++i) {
    species.Add(cloneTuneSpecies(*other.speciesAt(i)));
  }

  detachFromDirectories();
}

Tune& Tune::operator=(const Tune& other)
{
  if (this == &other) {
    return *this;
  }

  TObject::operator=(other);

  formatVersion = other.formatVersion;
  finalized = other.finalized;

  teacher = other.teacher;
  pythiaCard = other.pythiaCard;
  pythiaCardContent = other.pythiaCardContent;

  beamIdA = other.beamIdA;
  beamIdB = other.beamIdB;
  beamFrameType = other.beamFrameType;
  sqrtSNN = other.sqrtSNN;
  azimuthModel = other.azimuthModel;
  finalStatus = other.finalStatus;

  nEvents = other.nEvents;
  nAttempts = other.nAttempts;
  generationTimeSeconds = other.generationTimeSeconds;

  activityOverflowEvents = other.activityOverflowEvents;
  selectedMultiplicityOverflowEvents = other.selectedMultiplicityOverflowEvents;
  ptOverflowParticles = other.ptOverflowParticles;
  speciesMultiplicityOverflowEvents = other.speciesMultiplicityOverflowEvents;

  centralChargedCoverageMismatchEvents = other.centralChargedCoverageMismatchEvents;
  centralChargedCoverageMissingParticles = other.centralChargedCoverageMissingParticles;

  activityEtaMax = other.activityEtaMax;
  particleEtaMax = other.particleEtaMax;
  ptMax = other.ptMax;

  activityEdges = other.activityEdges;
  selectedMultiplicityEdges = other.selectedMultiplicityEdges;

  hNch = other.hNch;
  hNSelected = other.hNSelected;
  hNchSelected = other.hNchSelected;

  reconstructSparseLike(other.hNSelectedVsNch, hNSelectedVsNch);
  reconstructSparseLike(other.hNchSelectedVsNch, hNchSelectedVsNch);
  hEventsVsActivitySelected = other.hEventsVsActivitySelected;

  pNch = other.pNch;
  pNSelected = other.pNSelected;
  pNchSelected = other.pNchSelected;

  reconstructSparseLike(other.pNSelectedGivenNch, pNSelectedGivenNch);
  reconstructSparseLike(other.pNchSelectedGivenNch, pNchSelectedGivenNch);

  compositionTemplateCapPerPair = other.compositionTemplateCapPerPair;
  compositionPairNch = other.compositionPairNch;
  compositionPairNSelected = other.compositionPairNSelected;
  compositionPairOffsets = other.compositionPairOffsets;
  compositionPairEventsSeen = other.compositionPairEventsSeen;
  compositionCentralCounts = other.compositionCentralCounts;
  compositionOtherCounts = other.compositionOtherCounts;

  species.Delete();
  species.SetOwner(kTRUE);

  for (int i = 0; i < other.numberOfSpecies(); ++i) {
    species.Add(cloneTuneSpecies(*other.speciesAt(i)));
  }

  detachFromDirectories();

  return *this;
}

Tune::~Tune()
{
  species.Delete();
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
  formatVersion = kFormatVersion;
  finalized = false;

  teacher = "PYTHIA8";
  pythiaCard.clear();
  pythiaCardContent.clear();

  beamIdA = 0;
  beamIdB = 0;
  beamFrameType = 0;
  sqrtSNN = 0.0;
  azimuthModel = "uniform";
  finalStatus = 1;

  activityEtaMax = activityEtaMaxIn;
  particleEtaMax = particleEtaMaxIn;
  ptMax = ptMaxIn;
  activityEdges = activityEdgesIn;
  selectedMultiplicityEdges = selectedMultiplicityEdgesIn;

  hNch = TH1D("hNch", ";N_{ch};events", maxNch + 1, -0.5, maxNch + 0.5);

  hNSelected = TH1D("hNSelected", ";N_{selected};events", maxSelectedMultiplicity + 1, -0.5, maxSelectedMultiplicity + 0.5);

  hNchSelected = TH1D("hNchSelected", ";N_{ch}^{selected};events", maxNch + 1, -0.5, maxNch + 0.5);

  reconstructSparse(hNSelectedVsNch,
                    "hNSelectedVsNch",
                    ";N_{ch};N_{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxSelectedMultiplicity + 1,
                    -0.5,
                    maxSelectedMultiplicity + 0.5);

  reconstructSparse(hNchSelectedVsNch,
                    "hNchSelectedVsNch",
                    ";N_{ch};N_{ch}^{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5);

  const int nActivityBins = static_cast<int>(activityEdges.size()) - 1;
  const int nSelectedBins = static_cast<int>(selectedMultiplicityEdges.size()) - 1;

  hEventsVsActivitySelected = TH2D("hEventsVsActivitySelected", ";N_{ch};N_{selected}", nActivityBins, activityEdges.data(), nSelectedBins, selectedMultiplicityEdges.data());

  hNch.Sumw2();
  hNSelected.Sumw2();
  hNchSelected.Sumw2();

  pNch = TH1D();
  pNSelected = TH1D();
  pNchSelected = TH1D();

  reconstructSparse(pNSelectedGivenNch,
                    "pNSelectedGivenNch",
                    ";N_{ch};N_{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxSelectedMultiplicity + 1,
                    -0.5,
                    maxSelectedMultiplicity + 0.5);

  reconstructSparse(pNchSelectedGivenNch,
                    "pNchSelectedGivenNch",
                    ";N_{ch};N_{ch}^{selected}",
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5,
                    maxNch + 1,
                    -0.5,
                    maxNch + 0.5);

  compositionTemplateCapPerPair = 0;
  compositionPairNch.clear();
  compositionPairNSelected.clear();
  compositionPairOffsets.clear();
  compositionPairEventsSeen.clear();
  compositionCentralCounts.clear();
  compositionOtherCounts.clear();

  species.Delete();
  species.SetOwner(kTRUE);

  for (const int pdg : speciesPdg) {
    species.Add(new TuneSpecies(pdg,
                                nActivityBins,
                                activityEdges.data(),
                                nSelectedBins,
                                selectedMultiplicityEdges.data(),
                                maxSpeciesMultiplicity,
                                maxNch,
                                nPtBins,
                                ptMax,
                                nEtaBins,
                                particleEtaMax));
  }

  detachFromDirectories();
}

void Tune::finalize()
{
  pNch = hNch;
  pNSelected = hNSelected;
  pNchSelected = hNchSelected;

  pNch.SetName("pNch");
  pNSelected.SetName("pNSelected");
  pNchSelected.SetName("pNchSelected");

  normalize1D(pNch);
  normalize1D(pNSelected);
  normalize1D(pNchSelected);

  copyNormalizedSparseYSlices(hNSelectedVsNch, pNSelectedGivenNch);
  copyNormalizedSparseYSlices(hNchSelectedVsNch, pNchSelectedGivenNch);

  for (int i = 0; i < numberOfSpecies(); ++i) {
    speciesAt(i)->finalize();
  }

  normalizeSpeciesFractions(species, true);
  normalizeSpeciesFractions(species, false);

  finalized = true;

  detachFromDirectories();
}

void Tune::validate(bool requireFinalized) const
{
  if (formatVersion != kFormatVersion) {
    throw std::runtime_error("Ditto::Tune: unsupported tune format version");
  }

  if (teacher.empty() ||
      pythiaCardContent.empty()) {
    throw std::runtime_error("Ditto::Tune: missing teacher generator-card metadata");
  }

  if (beamIdA == 0 || beamIdB == 0) {
    throw std::runtime_error("Ditto::Tune: invalid beam IDs");
  }

  if (beamFrameType != 1) {
    throw std::runtime_error("Ditto::Tune: only Beams:frameType = 1 is supported");
  }

  if (!std::isfinite(sqrtSNN) || sqrtSNN <= 0.0) {
    throw std::runtime_error("Ditto::Tune: invalid sqrtSNN");
  }

  if (azimuthModel != "uniform") {
    throw std::runtime_error("Ditto::Tune: unsupported azimuth model");
  }

  if (finalStatus <= 0) {
    throw std::runtime_error("Ditto::Tune: finalStatus must be positive");
  }

  if (activityEtaMax <= 0.0) {
    throw std::runtime_error("Ditto::Tune: invalid activityEtaMax");
  }

  if (particleEtaMax <= 0.0) {
    throw std::runtime_error("Ditto::Tune: invalid particleEtaMax");
  }

  if (activityEdges.size() < 2) {
    throw std::runtime_error("Ditto::Tune: invalid activity edges");
  }

  if (selectedMultiplicityEdges.size() < 2) {
    throw std::runtime_error("Ditto::Tune: invalid selected-multiplicity edges");
  }

  for (std::size_t i = 1; i < activityEdges.size(); ++i) {
    if (activityEdges[i] <= activityEdges[i - 1]) {
      throw std::runtime_error("Ditto::Tune: activity edges are not strictly increasing");
    }
  }

  for (std::size_t i = 1; i < selectedMultiplicityEdges.size(); ++i) {
    if (selectedMultiplicityEdges[i] <= selectedMultiplicityEdges[i - 1]) {
      throw std::runtime_error("Ditto::Tune: selected-multiplicity edges are not strictly increasing");
    }
  }

  if (numberOfSpecies() <= 0) {
    throw std::runtime_error("Ditto::Tune: no species stored");
  }

  if (requireFinalized && !finalized) {
    throw std::runtime_error("Ditto::Tune: tune has not been finalized");
  }

  const int nActivityBins = static_cast<int>(activityEdges.size()) - 1;
  const int nSelectedBins = static_cast<int>(selectedMultiplicityEdges.size()) - 1;

  if (hEventsVsActivitySelected.GetNbinsX() != nActivityBins ||
      hEventsVsActivitySelected.GetNbinsY() != nSelectedBins) {
    throw std::runtime_error("Ditto::Tune: inconsistent event-composition binning");
  }

  if (requireFinalized) {
    if (pNSelectedGivenNch.GetNdimensions() != 2 ||
        pNchSelectedGivenNch.GetNdimensions() != 2 ||
        pNSelectedGivenNch.GetAxis(0)->GetNbins() != hNch.GetNbinsX() ||
        pNSelectedGivenNch.GetAxis(1)->GetNbins() != hNSelected.GetNbinsX() ||
        pNchSelectedGivenNch.GetAxis(0)->GetNbins() != hNch.GetNbinsX() ||
        pNchSelectedGivenNch.GetAxis(1)->GetNbins() != hNchSelected.GetNbinsX()) {
      throw std::runtime_error("Ditto::Tune: inconsistent sparse exact-multiplicity binning");
    }

    if (compositionPairNch.size() != compositionPairNSelected.size() ||
        compositionPairNch.size() != compositionPairEventsSeen.size()) {
      throw std::runtime_error("Ditto::Tune: inconsistent composition-pair metadata");
    }

    if (compositionPairOffsets.size() != compositionPairNch.size() + 1 ||
        compositionPairOffsets.empty() ||
        compositionPairOffsets.front() != 0) {
      throw std::runtime_error("Ditto::Tune: invalid composition-template offsets");
    }

    for (std::size_t i = 1; i < compositionPairOffsets.size(); ++i) {
      if (compositionPairOffsets[i] < compositionPairOffsets[i - 1]) {
        throw std::runtime_error("Ditto::Tune: composition-template offsets are not monotonic");
      }
    }

    for (std::size_t i = 1; i < compositionPairNch.size(); ++i) {
      if (compositionPairNch[i] < compositionPairNch[i - 1] ||
          (compositionPairNch[i] == compositionPairNch[i - 1] &&
           compositionPairNSelected[i] <= compositionPairNSelected[i - 1])) {
        throw std::runtime_error("Ditto::Tune: composition pairs are not strictly ordered");
      }
    }

    if (compositionTemplateCapPerPair == 0) {
      throw std::runtime_error("Ditto::Tune: invalid composition-template cap");
    }

    const std::uint64_t nTemplates =
      compositionPairOffsets.back();

    const std::uint64_t expectedCounts =
      nTemplates *
      static_cast<std::uint64_t>(numberOfSpecies());

    if (compositionCentralCounts.size() != expectedCounts ||
        compositionOtherCounts.size() != expectedCounts) {
      throw std::runtime_error("Ditto::Tune: inconsistent flattened composition-template size");
    }

    const std::size_t nSpecies =
      static_cast<std::size_t>(numberOfSpecies());

    for (std::size_t iPair = 0;
         iPair < compositionPairNch.size();
         ++iPair) {
      const std::uint64_t first =
        compositionPairOffsets[iPair];

      const std::uint64_t last =
        compositionPairOffsets[iPair + 1];

      const std::uint64_t stored =
        last - first;

      if (stored == 0 ||
          stored > compositionTemplateCapPerPair ||
          stored > compositionPairEventsSeen[iPair]) {
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
              compositionCentralCounts[static_cast<std::size_t>(base) + iSpecies]);

          const int other =
            static_cast<int>(
              compositionOtherCounts[static_cast<std::size_t>(base) + iSpecies]);

          centralTotal += central;
          selectedTotal += central + other;
        }

        if (centralTotal > compositionPairNch[iPair] ||
            selectedTotal != compositionPairNSelected[iPair]) {
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

    if (!std::isfinite(entry->mass) ||
        entry->mass < 0.0) {
      throw std::runtime_error("Ditto::Tune: invalid stored species mass");
    }

    if (entry->hCountVsActivity.GetNbinsX() != nActivityBins ||
        entry->hPtVsActivity.GetNbinsX() != nActivityBins ||
        entry->hEtaVsActivity.GetNbinsX() != nActivityBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent activity binning");
    }

    const int nExactNchBins =
      hNch.GetNbinsX();

    if (entry->hCentralChargedPtSumVsNch.GetNbinsX() != nExactNchBins ||
        entry->hCentralChargedPtCountVsNch.GetNbinsX() != nExactNchBins ||
        entry->hOtherPtSumVsNch.GetNbinsX() != nExactNchBins ||
        entry->hOtherPtCountVsNch.GetNbinsX() != nExactNchBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent exact-Nch pT-moment binning");
    }

    if (requireFinalized &&
        (entry->hCentralChargedMeanPtVsNch.GetNbinsX() != nExactNchBins ||
         entry->hOtherMeanPtVsNch.GetNbinsX() != nExactNchBins)) {
      throw std::runtime_error("Ditto::Tune: incomplete exact-Nch mean-pT tables");
    }

    if (entry->hCentralChargedCountVsActivitySelected.GetNbinsX() != nActivityBins ||
        entry->hCentralChargedCountVsActivitySelected.GetNbinsY() != nSelectedBins ||
        entry->hOtherCountVsActivitySelected.GetNbinsX() != nActivityBins ||
        entry->hOtherCountVsActivitySelected.GetNbinsY() != nSelectedBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent constrained-composition binning");
    }

    const int nSpeciesMultiplicityBins =
      entry->hCountVsActivity.GetNbinsY();

    if (entry->hCentralChargedMultiplicityVsActivity.GetNdimensions() != 2 ||
        entry->hOtherMultiplicityVsActivity.GetNdimensions() != 2 ||
        entry->hCentralChargedMultiplicityVsActivity.GetAxis(0)->GetNbins() != nActivityBins ||
        entry->hOtherMultiplicityVsActivity.GetAxis(0)->GetNbins() != nActivityBins ||
        entry->hCentralChargedMultiplicityVsActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins ||
        entry->hOtherMultiplicityVsActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins) {
      throw std::runtime_error("Ditto::Tune: inconsistent sparse multiplicity-diagnostic binning");
    }

    if (requireFinalized &&
        (entry->pCountGivenActivity.GetNbinsX() != nActivityBins ||
         entry->pPtGivenActivity.GetNbinsX() != nActivityBins ||
         entry->pEtaGivenActivity.GetNbinsX() != nActivityBins ||
         entry->pCentralChargedSpeciesGivenActivitySelected.GetNbinsX() != nActivityBins ||
         entry->pCentralChargedSpeciesGivenActivitySelected.GetNbinsY() != nSelectedBins ||
         entry->pOtherSpeciesGivenActivitySelected.GetNbinsX() != nActivityBins ||
         entry->pOtherSpeciesGivenActivitySelected.GetNbinsY() != nSelectedBins ||
         entry->pCentralChargedMultiplicityGivenActivity.GetNdimensions() != 2 ||
         entry->pOtherMultiplicityGivenActivity.GetNdimensions() != 2 ||
         entry->pCentralChargedMultiplicityGivenActivity.GetAxis(0)->GetNbins() != nActivityBins ||
         entry->pOtherMultiplicityGivenActivity.GetAxis(0)->GetNbins() != nActivityBins ||
         entry->pCentralChargedMultiplicityGivenActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins ||
         entry->pOtherMultiplicityGivenActivity.GetAxis(1)->GetNbins() != nSpeciesMultiplicityBins)) {
      throw std::runtime_error("Ditto::Tune: incomplete probability tables");
    }
  }
}

int Tune::activityClass(double nch) const
{
  if (activityEdges.size() < 2 ||
      nch < activityEdges.front() ||
      nch >= activityEdges.back()) {
    return -1;
  }

  const auto upper = std::upper_bound(activityEdges.begin(), activityEdges.end(), nch);

  return static_cast<int>(std::distance(activityEdges.begin(), upper) - 1);
}

int Tune::selectedMultiplicityClass(double nSelected) const
{
  if (selectedMultiplicityEdges.size() < 2 ||
      nSelected < selectedMultiplicityEdges.front() ||
      nSelected >= selectedMultiplicityEdges.back()) {
    return -1;
  }

  const auto upper = std::upper_bound(selectedMultiplicityEdges.begin(),
                                      selectedMultiplicityEdges.end(),
                                      nSelected);

  return static_cast<int>(std::distance(selectedMultiplicityEdges.begin(), upper) - 1);
}

int Tune::numberOfSpecies() const
{
  return species.GetEntriesFast();
}

int Tune::numberOfCompositionPairs() const
{
  return static_cast<int>(compositionPairNch.size());
}

std::uint64_t Tune::numberOfCompositionTemplates() const
{
  return compositionPairOffsets.empty()
           ? 0
           : compositionPairOffsets.back();
}

int Tune::compositionPairIndex(int nch, int nSelected) const
{
  const auto lower =
    std::lower_bound(
      compositionPairNch.begin(),
      compositionPairNch.end(),
      nch);

  std::size_t index =
    static_cast<std::size_t>(
      std::distance(compositionPairNch.begin(), lower));

  while (index < compositionPairNch.size() &&
         compositionPairNch[index] == nch) {
    if (compositionPairNSelected[index] == nSelected) {
      return static_cast<int>(index);
    }

    if (compositionPairNSelected[index] > nSelected) {
      break;
    }

    ++index;
  }

  return -1;
}

TuneSpecies* Tune::speciesAt(int index)
{
  return static_cast<TuneSpecies*>(species.UncheckedAt(index));
}

const TuneSpecies* Tune::speciesAt(int index) const
{
  return static_cast<const TuneSpecies*>(species.UncheckedAt(index));
}

TuneSpecies* Tune::findSpecies(int pdg)
{
  for (int i = 0; i < numberOfSpecies(); ++i) {
    auto* entry = speciesAt(i);

    if (entry && entry->pdg == pdg) {
      return entry;
    }
  }

  return nullptr;
}

const TuneSpecies* Tune::findSpecies(int pdg) const
{
  for (int i = 0; i < numberOfSpecies(); ++i) {
    const auto* entry = speciesAt(i);

    if (entry && entry->pdg == pdg) {
      return entry;
    }
  }

  return nullptr;
}

void Tune::detachFromDirectories()
{
  hNch.SetDirectory(nullptr);
  hNSelected.SetDirectory(nullptr);
  hNchSelected.SetDirectory(nullptr);

  hEventsVsActivitySelected.SetDirectory(nullptr);

  pNch.SetDirectory(nullptr);
  pNSelected.SetDirectory(nullptr);
  pNchSelected.SetDirectory(nullptr);


  species.SetOwner(kTRUE);

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
