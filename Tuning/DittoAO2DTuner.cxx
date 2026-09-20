///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoAO2DTuner.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/20
/// \brief  Standalone ROOT/AO2D source adapter for Ditto tune accumulation.
///

#include "DittoAO2DTuner.h"

#include "DittoTune.h"

#include <TClass.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGrid.h>
#include <TKey.h>
#include <TList.h>
#include <TParticlePDG.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Ditto
{
namespace
{

constexpr std::uint8_t kProducedByTransport = 0x1;
constexpr std::uint8_t kPhysicalPrimary = 0x4;
constexpr float kMaxPhysicalPrimaryRadius = 5.f; // cm, O2 definition

bool startsWith(const std::string& value, const std::string& prefix)
{
  return value.size() >= prefix.size() && value.compare(0, prefix.size(), prefix) == 0;
}

/// Extract the HepMC status from the O2 MCGenStatusEncoding without requiring
/// AliceO2 headers. The encoding reserves the upper three bits as marker 5 and
/// the low nine bits for the signed HepMC status. Old AO2Ds simply store the
/// HepMC status directly and therefore pass through unchanged.
int hepMCStatusCode(int statusCode)
{
  const std::uint32_t bits = static_cast<std::uint32_t>(statusCode);
  constexpr std::uint32_t encodedMarker = 5u;
  constexpr std::uint32_t markerShift = 29u;

  if ((bits >> markerShift) != encodedMarker) {
    return statusCode;
  }

  constexpr std::uint32_t hepMask = 0x1ffu; // signed 9-bit value
  std::uint32_t raw = bits & hepMask;
  if ((raw & 0x100u) != 0u) {
    raw |= ~hepMask; // sign-extend bit 8
  }
  return static_cast<int>(static_cast<std::int32_t>(raw));
}

std::uint8_t protectedFlags(std::uint8_t flags, float vx, float vy)
{
  if (std::hypot(vx, vy) > kMaxPhysicalPrimaryRadius && (flags & kPhysicalPrimary) != 0u) {
    flags = static_cast<std::uint8_t>(flags & ~kPhysicalPrimary);
  }
  return flags;
}

bool isSelected(AO2DParticleSelection selection,
                std::uint8_t flags,
                int statusCode,
                float vx,
                float vy)
{
  switch (selection) {
    case AO2DParticleSelection::GeneratorFinal:
      return (flags & kProducedByTransport) == 0u &&
             hepMCStatusCode(statusCode) == 1;

    case AO2DParticleSelection::PhysicalPrimary:
      return (protectedFlags(flags, vx, vy) & kPhysicalPrimary) != 0u;
  }

  return false;
}

/// Return the newest tree in a directory whose name starts with prefix.
/// This accepts both unversioned and versioned AO2D tables, e.g.
/// O2mcparticle, O2mcparticle_000, O2mcparticle_001.
TTree* findTree(TDirectory& directory, const std::string& prefix)
{
  std::string bestName;

  TIter next(directory.GetListOfKeys());
  while (auto* object = next()) {
    auto* key = dynamic_cast<TKey*>(object);
    if (!key) {
      continue;
    }

    const std::string name = key->GetName();
    if (!startsWith(name, prefix)) {
      continue;
    }

    TClass* cls = TClass::GetClass(key->GetClassName());
    if (!cls || !cls->InheritsFrom(TTree::Class())) {
      continue;
    }

    if (name == prefix) {
      bestName = name;
      break;
    }
    if (startsWith(name, prefix + "_") && name > bestName) {
      bestName = name;
    }
  }

  if (bestName.empty()) {
    throw std::runtime_error("Ditto::AO2DTuner: no AO2D table matching prefix " + prefix +
                             " found in " + directory.GetPath());
    return nullptr;
  }
  // Final check on the naming, ideally the found name should be equal to prefix or start with prefix + "_". This is a sanity check to avoid accidentally picking up an unrelated tree.
  if (bestName != prefix && !startsWith(bestName, prefix + "_")) {
    throw std::runtime_error("Ditto::AO2DTuner: found AO2D table " + bestName +
                             " in " + directory.GetPath() +
                             " does not match expected prefix " + prefix);
  }
  return dynamic_cast<TTree*>(directory.Get(bestName.c_str()));
}

std::vector<TDirectory*> dataframeDirectories(TFile& file)
{
  std::vector<std::pair<std::string, TDirectory*>> found;

  TIter next(file.GetListOfKeys());
  while (auto* object = next()) {
    auto* key = dynamic_cast<TKey*>(object);
    if (!key) {
      continue;
    }

    const std::string name = key->GetName();
    if (!startsWith(name, "DF_")) {
      continue;
    }

    TClass* cls = TClass::GetClass(key->GetClassName());
    if (!cls || !cls->InheritsFrom(TDirectory::Class())) {
      continue;
    }

    if (auto* dir = file.GetDirectory(name.c_str())) {
      found.emplace_back(name, dir);
    }
  }

  std::sort(found.begin(), found.end(), [](const auto& a, const auto& b) {
    return a.first < b.first;
  });

  std::vector<TDirectory*> result;
  result.reserve(found.size());
  for (const auto& [name, dir] : found) {
    (void)name;
    result.push_back(dir);
  }

  // Some small/private AO2Ds may store the tables at file root. Support that
  // layout as a convenience if there are no DF_* directories.
  if (result.empty() && findTree(file, "O2mcparticle") != nullptr) {
    result.push_back(&file);
  }

  return result;
}

double etaFromMomentum(float px, float py, float pz)
{
  const double pt = std::hypot(static_cast<double>(px), static_cast<double>(py));
  if (pt > 0.0) {
    return std::asinh(static_cast<double>(pz) / pt);
  }
  if (pz > 0.f) {
    return 100.0;
  }
  if (pz < 0.f) {
    return -100.0;
  }
  return 0.0;
}

bool isNuclearPdg(int pdg)
{
  const std::int64_t absPdg = std::llabs(static_cast<long long>(pdg));
  return absPdg >= 1000000000LL;
}

int nuclearChargeType(int pdg)
{
  // PDG ion convention: +/-10LZZZAAAI.
  const std::int64_t absPdg = std::llabs(static_cast<long long>(pdg));
  const int z = static_cast<int>((absPdg / 10000LL) % 1000LL);
  const int sign = pdg < 0 ? -1 : 1;
  return 3 * sign * z;
}

class ParticleDatabase
{
 public:
  int chargeType(int pdg)
  {
    const auto found = mChargeCache.find(pdg);
    if (found != mChargeCache.end()) {
      return found->second;
    }

    int charge = 0;
    if (const auto* particle = TDatabasePDG::Instance()->GetParticle(pdg)) {
      charge = static_cast<int>(std::lround(particle->Charge()));
    } else if (isNuclearPdg(pdg)) {
      charge = nuclearChargeType(pdg);
    } else {
      throw std::runtime_error("unknown PDG code " + std::to_string(pdg));
    }

    mChargeCache.emplace(pdg, charge);
    return charge;
  }

  ParticleDefinition definition(int pdg)
  {
    if (const auto* particle = TDatabasePDG::Instance()->GetParticle(pdg)) {
      return ParticleDefinition{pdg,
                                particle->GetName(),
                                particle->Mass(),
                                static_cast<int>(std::lround(particle->Charge()))};
    }

    if (isNuclearPdg(pdg)) {
      throw std::runtime_error(
        "configured nuclear species " + std::to_string(pdg) +
        " is not present in ROOT TDatabasePDG; provide a ROOT PDG entry before constructing AO2DTuner");
    }

    throw std::runtime_error("configured species has unknown PDG code " +
                             std::to_string(pdg));
  }

 private:
  std::unordered_map<int, int> mChargeCache;
};

} // namespace

struct AO2DTunerImpl {
  explicit AO2DTunerImpl(const AO2DTunerConfig& cfg)
    : config(cfg)
  {
    validateConfig();
    configureAccumulator();
  }

  AO2DTunerConfig config;
  ParticleDatabase particleDatabase;
  std::unique_ptr<TuneAccumulator> accumulator;
  InputEvent eventBuffer;

  std::uint64_t processedParticles = 0;
  std::uint64_t selectedParticles = 0;
  std::uint64_t unknownPdgParticles = 0;
  bool ran = false;

  void validateConfig() const
  {
    if (config.mPythiaCard.empty()) {
      throw std::invalid_argument("Ditto::AO2DTuner: mPythiaCard must not be empty");
    }
    // Check that the card points at an existing file
    if (!std::ifstream(config.mPythiaCard).good()) {
      throw std::invalid_argument("Ditto::AO2DTuner: mPythiaCard file does not exist or is not readable: " + config.mPythiaCard);
    }
    if (config.mInputFiles.empty()) {
      throw std::invalid_argument("Ditto::AO2DTuner: mInputFiles must not be empty");
    }
  }

  void configureAccumulator()
  {
    std::vector<ParticleDefinition> definitions;
    definitions.reserve(config.mSpecies.size());
    for (const int pdg : config.mSpecies) {
      definitions.push_back(particleDatabase.definition(pdg));
    }

    accumulator = std::make_unique<TuneAccumulator>(config, definitions);

    auto& outputTune = accumulator->tune();
    outputTune.mTeacher = "AO2D-MC";
    outputTune.importPythiaCard(config.mPythiaCard);
  }

  bool reachedEventLimit() const
  {
    return config.mMaxEvents > 0 && accumulator->processedEvents() >= config.mMaxEvents;
  }

  void submitEvent()
  {
    if (reachedEventLimit()) {
      return;
    }

    accumulator->processEvent(eventBuffer);
    eventBuffer.clear();

    const auto n = accumulator->processedEvents();
    if (config.mProgressEvery > 0 && n % config.mProgressEvery == 0) {
      std::cout << "Ditto AO2D tuner: " << n << " events processed\n";
    }
  }

  bool addParticle(int pdg,
                   int statusCode,
                   std::uint8_t flags,
                   float px,
                   float py,
                   float pz,
                   float vx,
                   float vy)
  {
    ++processedParticles;

    if (!isSelected(config.mParticleSelection, flags, statusCode, vx, vy)) {
      return false;
    }

    int charge = 0;
    try {
      charge = particleDatabase.chargeType(pdg);
    } catch (const std::exception&) {
      ++unknownPdgParticles;
      if (config.mIgnoreUnknownPdg) {
        return false;
      }
      throw;
    }

    const double pt = std::hypot(static_cast<double>(px), static_cast<double>(py));
    const double eta = etaFromMomentum(px, py, pz);

    eventBuffer.push_back(InputParticle{pdg, pt, eta, charge});
    ++selectedParticles;
    return true;
  }

  void processDirectory(TDirectory& directory, const std::string& fileName)
  {
    auto* particleTree = findTree(directory, "O2mcparticle");
    auto* collisionTree = findTree(directory, "O2mccollision");

    if (!particleTree || !collisionTree) {
      throw std::runtime_error("Ditto::AO2DTuner: missing O2mcparticle/O2mccollision table in " + fileName + ":" + directory.GetPath());
    }

    const Long64_t nCollisions = collisionTree->GetEntries();
    if (nCollisions < 0) {
      throw std::runtime_error("Ditto::AO2DTuner: invalid MC collision count");
    }

    TTreeReader reader(particleTree);
    TTreeReaderValue<int> mcCollisionId(reader, "fIndexMcCollisions");
    TTreeReaderValue<int> pdgCode(reader, "fPdgCode");
    TTreeReaderValue<int> statusCode(reader, "fStatusCode");
    TTreeReaderValue<UChar_t> flags(reader, "fFlags");
    TTreeReaderValue<float> px(reader, "fPx");
    TTreeReaderValue<float> py(reader, "fPy");
    TTreeReaderValue<float> pz(reader, "fPz");
    TTreeReaderValue<float> vx(reader, "fVx");
    TTreeReaderValue<float> vy(reader, "fVy");

    Long64_t currentCollision = 0;
    eventBuffer.clear();

    while (reader.Next()) {
      if (reachedEventLimit()) {
        return;
      }

      const int collision = *mcCollisionId;
      if (collision < 0 || collision >= nCollisions) {
        throw std::runtime_error("Ditto::AO2DTuner: MC particle refers to invalid collision index " + std::to_string(collision) +
                                 " while nCollisions = " + std::to_string(nCollisions));
      }
      if (collision < currentCollision) {
        throw std::runtime_error(
          "Ditto::AO2DTuner: O2mcparticle is not ordered by MC collision index; "
          "streaming grouping cannot be used safely");
      }

      while (currentCollision < collision) {
        submitEvent();
        if (reachedEventLimit()) {
          return;
        }
        ++currentCollision;
      }

      addParticle(*pdgCode,
                  *statusCode,
                  static_cast<std::uint8_t>(*flags),
                  *px,
                  *py,
                  *pz,
                  *vx,
                  *vy);
    }

    while (currentCollision < nCollisions && !reachedEventLimit()) {
      submitEvent();
      ++currentCollision;
    }
  }

  void processFile(const std::string& fileName)
  {
    if (startsWith(fileName, "alien://") && !gGrid) {
      if (TGrid::Connect("alien://") == nullptr) {
        throw std::runtime_error("Ditto::AO2DTuner: could not connect to AliEn");
      }
    }
    std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
    if (!file || file->IsZombie()) {
      throw std::runtime_error("Ditto::AO2DTuner: could not open AO2D file: " + fileName);
    }

    const auto directories = dataframeDirectories(*file);
    if (directories.empty()) {
      throw std::runtime_error("Ditto::AO2DTuner: no DF_* directories containing O2mcparticle found in " + fileName);
    }

    for (auto* directory : directories) {
      if (reachedEventLimit()) {
        break;
      }
      processDirectory(*directory, fileName);
    }
  }

  void run()
  {
    if (ran) {
      throw std::runtime_error("Ditto::AO2DTuner: run() may only be called once");
    }
    if (config.mInputFiles.empty()) {
      throw std::invalid_argument("Ditto::AO2DTuner: no AO2D input files configured");
    }

    ran = true;
    const auto start = std::chrono::steady_clock::now();

    int counter = 0;
    for (const auto& fileName : config.mInputFiles) {
      if (reachedEventLimit()) {
        break;
      }
      const auto timeBefore = std::chrono::steady_clock::now();
      std::cout << "Ditto AO2D tuner: reading file " << counter + 1 << "/" << config.mInputFiles.size() << ": " << fileName << "\n";
      processFile(fileName);
      const auto timeAfter = std::chrono::steady_clock::now();
      const auto singleFileElapsed = std::chrono::duration<double>(timeAfter - timeBefore).count();
      std::cout << "\tDitto AO2D tuner: finished reading file " << counter + 1 << "/" << config.mInputFiles.size() << " in " << singleFileElapsed << "s\n";
      const auto totalElapsed = std::chrono::duration<double>(timeAfter - start).count();
      const auto totalEvents = accumulator->processedEvents();
      const auto secondsPerEvent = totalEvents > 0 ? totalElapsed / static_cast<double>(totalEvents) : 0.0;
      const auto secondsPerFile = (counter + 1) > 0 ? totalElapsed / static_cast<double>(counter + 1) : 0.0;
      std::cout << "\t\tEvent rate: " << 1.0 / secondsPerEvent << " events/s.";
      std::cout << "\tFile rate: " << 1.0 / secondsPerFile << " files/s.";
      if (config.mMaxEvents > 0) {
        std::cout << "\tETA: " << (config.mMaxEvents - totalEvents) * secondsPerEvent << " s\n";
      } else {
        std::cout << "\tETA: " << (config.mInputFiles.size() - counter - 1) * secondsPerFile << " s\n";
      }
      counter++;
    }

    accumulator->finalize();

    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    std::cout << "\nDitto AO2D tuning complete\n"
              << "  events             : " << accumulator->processedEvents() << "\n"
              << "  particles read     : " << processedParticles << "\n"
              << "  particles selected : " << selectedParticles << "\n"
              << "  unknown PDGs       : " << unknownPdgParticles << "\n"
              << "  wall time          : " << elapsed << " s\n"
              << "  composition pairs  : " << accumulator->tune().numberOfCompositionPairs() << "\n"
              << "  templates stored   : " << accumulator->tune().numberOfCompositionTemplates() << "\n";

    if (accumulator->activityOverflowEvents() > 0) {
      std::cout << "  WARNING activity overflow events: " << accumulator->activityOverflowEvents() << "\n";
    }
    if (accumulator->ptOverflowParticles() > 0) {
      std::cout << "  WARNING pT overflow particles: " << accumulator->ptOverflowParticles() << "\n";
    }
    if (accumulator->speciesMultiplicityOverflowEvents() > 0) {
      std::cout << "  WARNING species-count overflow fills: " << accumulator->speciesMultiplicityOverflowEvents() << "\n";
    }
  }
};

AO2DTuner::AO2DTuner(const AO2DTunerConfig& config)
  : mImpl(new AO2DTunerImpl(config))
{
}

AO2DTuner::~AO2DTuner()
{
  delete mImpl;
  mImpl = nullptr;
}

void AO2DTuner::addFile(const std::string& fileName)
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::AO2DTuner: invalid implementation");
  }
  if (mImpl->ran) {
    throw std::runtime_error("Ditto::AO2DTuner: cannot add files after run()");
  }
  mImpl->config.mInputFiles.push_back(fileName);
}

void AO2DTuner::run()
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::AO2DTuner: invalid implementation");
  }
  mImpl->run();
}

void AO2DTuner::save(const std::string& fileName) const
{
  if (!mImpl || !mImpl->ran || !mImpl->accumulator->finalized()) {
    throw std::runtime_error("Ditto::AO2DTuner: run() must finish before save()");
  }
  mImpl->accumulator->tune().save(fileName);
}

const AO2DTunerConfig& AO2DTuner::config() const
{
  if (!mImpl) {
    throw std::runtime_error("Ditto::AO2DTuner: invalid implementation");
  }
  return mImpl->config;
}

std::uint64_t AO2DTuner::processedEvents() const
{
  return mImpl && mImpl->accumulator ? mImpl->accumulator->processedEvents() : 0;
}

std::uint64_t AO2DTuner::processedParticles() const
{
  return mImpl ? mImpl->processedParticles : 0;
}

std::uint64_t AO2DTuner::selectedParticles() const
{
  return mImpl ? mImpl->selectedParticles : 0;
}

std::uint64_t AO2DTuner::unknownPdgParticles() const
{
  return mImpl ? mImpl->unknownPdgParticles : 0;
}

const Tune& AO2DTuner::tune() const
{
  if (!mImpl || !mImpl->accumulator) {
    throw std::runtime_error("Ditto::AO2DTuner: accumulator is not initialized");
  }
  return mImpl->accumulator->tune();
}

} // namespace Ditto
