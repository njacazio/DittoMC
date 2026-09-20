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

#include <TChain.h>
#include <TClass.h>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGrid.h>
#include <TKey.h>
#include <TParticlePDG.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
      return (flags & kProducedByTransport) == 0u && hepMCStatusCode(statusCode) == 1;
    case AO2DParticleSelection::PhysicalPrimary:
      return (protectedFlags(flags, vx, vy) & kPhysicalPrimary) != 0u;
  }
  return false;
}

/// Return the AO2D tree name matching prefix in one directory.
///
/// The unversioned table is preferred when present. Otherwise the highest
/// lexicographic version is used, e.g. O2mcparticle_001 over O2mcparticle_000.
std::string findTreeName(TDirectory& directory, const std::string& prefix)
{
  std::string bestName;

  TIter next(directory.GetListOfKeys());
  while (auto* object = next()) {
    auto* key = dynamic_cast<TKey*>(object);
    if (!key) {
      continue;
    }

    const std::string name = key->GetName();
    if (name != prefix && !startsWith(name, prefix + "_")) {
      continue;
    }

    TClass* cls = TClass::GetClass(key->GetClassName());
    if (!cls || !cls->InheritsFrom(TTree::Class())) {
      continue;
    }

    if (name == prefix) {
      return name;
    }

    if (name > bestName) {
      bestName = name;
    }
  }

  if (bestName.empty()) {
    throw std::runtime_error("Ditto::AO2DTuner: no AO2D table matching prefix " +
                             prefix + " found in " + directory.GetPath());
  }

  return bestName;
}

/// Return DF_* directory names in deterministic order.
///
/// An empty string is used as a sentinel for small/private AO2Ds that store
/// the tables directly at file root.
std::vector<std::string> dataframeDirectoryNames(TFile& file)
{
  std::vector<std::string> result;

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

    result.push_back(name);
  }

  std::sort(result.begin(), result.end());
  if (result.empty()) {
    result.emplace_back();
  }

  return result;
}

TDirectory& dataframeDirectory(TFile& file, const std::string& name)
{
  if (name.empty()) {
    return file;
  }

  auto* directory = file.GetDirectory(name.c_str());
  if (!directory) {
    throw std::runtime_error("Ditto::AO2DTuner: could not open dataframe directory " +
                             name + " in " + file.GetName());
  }
  return *directory;
}

TTree& exactTree(TDirectory& directory,
                 const std::string& treeName,
                 const std::string& fileName)
{
  auto* tree = dynamic_cast<TTree*>(directory.Get(treeName.c_str()));
  if (!tree) {
    throw std::runtime_error(
      "Ditto::AO2DTuner: expected AO2D table " + treeName + " in " + fileName +
      ":" + directory.GetPath() +
      ". All input AO2Ds must use the schema discovered from the first file.");
  }
  return *tree;
}

std::string treePath(const std::string& directoryName, const std::string& treeName)
{
  return directoryName.empty() ? treeName : directoryName + "/" + treeName;
}

void ensureGridConnection(const std::string& fileName)
{
  if (startsWith(fileName, "alien://") && !gGrid) {
    if (TGrid::Connect("alien://") == nullptr) {
      throw std::runtime_error("Ditto::AO2DTuner: could not connect to AliEn");
    }
  }
}

std::unique_ptr<TFile> openAO2D(const std::string& fileName)
{
  ensureGridConnection(fileName);

  std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Ditto::AO2DTuner: could not open AO2D file: " + fileName);
  }
  return file;
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
        " is not present in ROOT TDatabasePDG; provide a ROOT PDG entry before "
        "constructing AO2DTuner");
    }

    throw std::runtime_error("configured species has unknown PDG code " +
                             std::to_string(pdg));
  }

 private:
  std::unordered_map<int, int> mChargeCache;
};

} // namespace

struct AO2DTunerImpl {
  struct InputChunk {
    InputChunk(std::string file, std::string directory, Long64_t nParticles, Long64_t nCollisions)
      : fileName(std::move(file)),
        directoryName(std::move(directory)),
        nParticleEntries(nParticles),
        nCollisions(nCollisions),
        fileIndex(nextFileIndex++)
    {
    }
    const std::string fileName;
    const std::string directoryName;
    const Long64_t nParticleEntries = 0;
    const Long64_t nCollisions = 0;

   private:
    static std::size_t nextFileIndex;
    const std::size_t fileIndex = 0;
  };

  explicit AO2DTunerImpl(const AO2DTunerConfig& cfg) : config(cfg)
  {
    validateConfig();
    configureAccumulator();
  }

  AO2DTunerConfig config;
  ParticleDatabase particleDatabase;
  std::unique_ptr<TuneAccumulator> accumulator;
  std::unique_ptr<TChain> particleChain;
  std::vector<InputChunk> chunks;
  InputEvent eventBuffer;

  std::string particleTreeName;
  std::string collisionTreeName;

  Long64_t totalParticleEntries = 0;
  std::uint64_t totalInputEvents = 0;
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

  std::string formatDuration(double seconds)
  {
    if (!std::isfinite(seconds) || seconds < 0.0) {
      return "--:--:--";
    }

    const auto total = static_cast<std::uint64_t>(seconds);

    const auto hours = total / 3600;
    const auto minutes = (total % 3600) / 60;
    const auto secs = total % 60;

    char buffer[32];
    std::snprintf(
      buffer,
      sizeof(buffer),
      "%02llu:%02llu:%02llu",
      static_cast<unsigned long long>(hours),
      static_cast<unsigned long long>(minutes),
      static_cast<unsigned long long>(secs));

    return buffer;
  }

  void discoverSchema(TFile& file, const std::vector<std::string>& directories)
  {
    if (directories.empty()) {
      throw std::runtime_error("Ditto::AO2DTuner: no dataframe directories found in " + std::string(file.GetName()));
    }

    auto& firstDirectory = dataframeDirectory(file, directories.front());
    particleTreeName = findTreeName(firstDirectory, "O2mcparticle");
    collisionTreeName = findTreeName(firstDirectory, "O2mccollision");

    std::cout << "Ditto AO2D tuner: schema from first input\n"
              << "  particle table  : " << particleTreeName << "\n"
              << "  collision table : " << collisionTreeName << "\n";
  }

  void discoverInputSchema()
  {
    auto file = openAO2D(config.mInputFiles.front());
    const auto directories = dataframeDirectoryNames(*file);
    discoverSchema(*file, directories);
  }

  void appendFileToChain(const std::string& fileName, std::size_t fileIndex)
  {
    auto file = openAO2D(fileName);
    const auto directories = dataframeDirectoryNames(*file);

    for (const auto& directoryName : directories) {
      auto& directory = dataframeDirectory(*file, directoryName);
      auto& particleTree = exactTree(directory, particleTreeName, fileName);
      auto& collisionTree = exactTree(directory, collisionTreeName, fileName);

      const Long64_t nParticles = particleTree.GetEntries();
      const Long64_t nCollisions = collisionTree.GetEntries();
      if (nParticles < 0 || nCollisions < 0) {
        throw std::runtime_error("Ditto::AO2DTuner: invalid AO2D entry count in " + fileName + ":" + directory.GetPath());
      }

      chunks.push_back(InputChunk{fileName,
                                  directoryName,
                                  nParticles,
                                  nCollisions});

      if (static_cast<unsigned long long>(nCollisions) > std::numeric_limits<std::uint64_t>::max() - totalInputEvents) {
        throw std::overflow_error("Ditto::AO2DTuner: total MC collision count overflow");
      }
      totalInputEvents += static_cast<std::uint64_t>(nCollisions);

      // Empty particle tables still correspond to valid empty MC collisions,
      // represented by the chunk metadata. They do not need a TChain element.
      if (nParticles == 0) {
        continue;
      }

      const std::string path = treePath(directoryName, particleTreeName);
      particleChain->AddFile(fileName.c_str(), nParticles, path.c_str());
      totalParticleEntries += nParticles;
    }
  }

  void buildInputChain(std::size_t firstFile, std::size_t lastFile)
  {
    chunks.clear();

    particleChain = std::make_unique<TChain>(particleTreeName.c_str());

    totalParticleEntries = 0;
    totalInputEvents = 0;

    for (std::size_t i = firstFile; i < lastFile; ++i) {
      std::cout << "Ditto AO2D tuner: indexing file "
                << i + 1 << "/" << config.mInputFiles.size()
                << ": " << config.mInputFiles[i] << "\n";

      appendFileToChain(config.mInputFiles[i], i);
    }

    if (!particleChain || chunks.empty()) {
      throw std::runtime_error("Ditto::AO2DTuner: no AO2D input chunks found");
    }

    std::cout << "Ditto AO2D tuner: batch contains "
              << chunks.size() << " dataframe chunks, "
              << totalParticleEntries << " MC-particle entries and "
              << totalInputEvents << " MC collisions\n";
  }

  void processEmptyInput()
  {
    for (const auto& chunk : chunks) {
      eventBuffer.clear();
      for (Long64_t collision = 0; collision < chunk.nCollisions && !reachedEventLimit(); ++collision) {
        submitEvent();
      }
      if (reachedEventLimit()) {
        return;
      }
    }
  }

  void processChain()
  {
    if (totalParticleEntries == 0) {
      processEmptyInput();
      return;
    }

    // The particle chain is the hot I/O path. Read only the branches used by
    // the tuner; this matters substantially for large AO2Ds and remote input.
    particleChain->SetBranchStatus("*", false);
    particleChain->SetBranchStatus("fIndexMcCollisions", true);
    particleChain->SetBranchStatus("fPdgCode", true);
    particleChain->SetBranchStatus("fStatusCode", true);
    particleChain->SetBranchStatus("fFlags", true);
    particleChain->SetBranchStatus("fPx", true);
    particleChain->SetBranchStatus("fPy", true);
    particleChain->SetBranchStatus("fPz", true);
    particleChain->SetBranchStatus("fVx", true);
    particleChain->SetBranchStatus("fVy", true);

    TTreeReader reader(particleChain.get());
    TTreeReaderValue<int> mcCollisionId(reader, "fIndexMcCollisions");
    TTreeReaderValue<int> pdgCode(reader, "fPdgCode");
    TTreeReaderValue<int> statusCode(reader, "fStatusCode");
    TTreeReaderValue<UChar_t> flags(reader, "fFlags");
    TTreeReaderValue<float> px(reader, "fPx");
    TTreeReaderValue<float> py(reader, "fPy");
    TTreeReaderValue<float> pz(reader, "fPz");
    TTreeReaderValue<float> vx(reader, "fVx");
    TTreeReaderValue<float> vy(reader, "fVy");

    Long64_t consumedParticleEntries = 0;

    for (const auto& chunk : chunks) {
      if (reachedEventLimit()) {
        return;
      }

      Long64_t currentCollision = 0;
      eventBuffer.clear();

      for (Long64_t localEntry = 0; localEntry < chunk.nParticleEntries; ++localEntry) {
        if (reachedEventLimit()) {
          return;
        }

        if (!reader.Next()) {
          throw std::runtime_error("Ditto::AO2DTuner: TChain ended before the indexed AO2D particle count was reached while reading " +
                                   chunk.fileName + ":" +
                                   (chunk.directoryName.empty() ? std::string("/") : chunk.directoryName));
        }
        ++consumedParticleEntries;

        const int collision = *mcCollisionId;
        if (collision < 0 || collision >= chunk.nCollisions) {
          throw std::runtime_error("Ditto::AO2DTuner: MC particle refers to invalid collision index " +
                                   std::to_string(collision) + " in " + chunk.fileName + ":" +
                                   (chunk.directoryName.empty() ? std::string("/") : chunk.directoryName) +
                                   " while nCollisions = " + std::to_string(chunk.nCollisions));
        }

        if (collision < currentCollision) {
          throw std::runtime_error("Ditto::AO2DTuner: " + particleTreeName +
                                   " is not ordered by MC collision index in " + chunk.fileName + ":" +
                                   (chunk.directoryName.empty() ? std::string("/") : chunk.directoryName) +
                                   "; streaming grouping cannot be used safely");
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

      // fIndexMcCollisions is local to each DF. Flush the remaining collisions
      // here and reset before the next chunk rather than treating the TChain as
      // one global collision-index space.
      while (currentCollision < chunk.nCollisions && !reachedEventLimit()) {
        submitEvent();
        ++currentCollision;
      }
    }

    if (!reachedEventLimit() && consumedParticleEntries != totalParticleEntries) {
      throw std::runtime_error("Ditto::AO2DTuner: particle-chain accounting mismatch: consumed " +
                               std::to_string(consumedParticleEntries) + " entries, expected " +
                               std::to_string(totalParticleEntries));
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

    //
    // Detect the AO2D schema once, from the first input file.
    //
    discoverInputSchema();

    const std::size_t nFiles = config.mInputFiles.size();
    const std::size_t batchSize = config.mFileBatchSize > 0 ? config.mFileBatchSize : nFiles;

    for (std::size_t first = 0; first < nFiles && !reachedEventLimit(); first += batchSize) {

      const std::size_t last = std::min(first + batchSize, nFiles);

      std::cout << "\nDitto AO2D tuner: processing file batch "
                << first + 1 << "-" << last
                << "/" << nFiles << "\n";

      buildInputChain(first, last);

      processChain();

      //
      // The chain and its file handles are no longer needed.
      //
      particleChain.reset();
      chunks.clear();

      //
      // Overall progress / ETA.
      //
      const std::size_t completedFiles = last;

      const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

      const double secondsPerFile = completedFiles > 0 ? elapsed / static_cast<double>(completedFiles) : 0.0;

      const double eta = secondsPerFile * static_cast<double>(nFiles - completedFiles);

      const double fraction = static_cast<double>(completedFiles) / static_cast<double>(nFiles);

      std::cout << "Ditto AO2D tuner: "
                << completedFiles << "/" << nFiles
                << " files complete"
                << " (" << 100.0 * fraction << "%)"
                << ", elapsed " << formatDuration(elapsed)
                << ", ETA " << formatDuration(eta)
                << "\n";
    }

    accumulator->finalize();

    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const auto nEvents = accumulator->processedEvents();
    const double eventRate = elapsed > 0.0 ? static_cast<double>(nEvents) / elapsed : 0.0;

    std::cout << "\nDitto AO2D tuning complete\n"
              << " events               : " << nEvents << "\n"
              << " particles read       : " << processedParticles << "\n"
              << " particles selected   : " << selectedParticles << "\n"
              << " unknown PDGs         : " << unknownPdgParticles << "\n"
              << " wall time            : " << elapsed << " s\n"
              << " event rate           : " << eventRate << " events/s\n"
              << " composition pairs    : "
              << accumulator->tune().numberOfCompositionPairs() << "\n"
              << " templates stored     : "
              << accumulator->tune().numberOfCompositionTemplates() << "\n";

    if (accumulator->activityOverflowEvents() > 0) {
      std::cout << " WARNING activity overflow events: " << accumulator->activityOverflowEvents() << "\n";
    }

    if (accumulator->ptOverflowParticles() > 0) {
      std::cout << " WARNING pT overflow particles: " << accumulator->ptOverflowParticles() << "\n";
    }

    if (accumulator->speciesMultiplicityOverflowEvents() > 0) {
      std::cout << " WARNING species-count overflow fills: " << accumulator->speciesMultiplicityOverflowEvents() << "\n";
    }
  }
};

std::size_t AO2DTunerImpl::InputChunk::nextFileIndex = 0;

AO2DTuner::AO2DTuner(const AO2DTunerConfig& config) : mImpl(new AO2DTunerImpl(config))
{
}

AO2DTuner::~AO2DTuner()
{
  delete mImpl;
  mImpl = nullptr;
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
