#include "DittoTune.h"
#include "DittoTuneAccumulator.h"

#include <TTree.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void require(bool condition, const std::string& message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void requireClose(double value, double expected, double tolerance, const std::string& message)
{
  if (std::abs(value - expected) > tolerance) {
    throw std::runtime_error(message + ": got " + std::to_string(value) + ", expected " + std::to_string(expected));
  }
}

} // namespace

int main()
{
  try {
    using namespace Ditto;

    TuneAccumulatorConfig cfg;

    cfg.mActivityEtaMax = 0.5;
    cfg.mParticleEtaMax = 2.0;

    cfg.mActivityEdges = {0.0, 5.0, 10.0};

    cfg.mSelectedMultiplicityEdges = {0.0, 2.0, 5.0, 10.0};

    cfg.mMaxNch = 20;
    cfg.mMaxSelectedMultiplicity = 20;
    cfg.mMaxSpeciesMultiplicity = 10;

    cfg.mNPtBins = 20;
    cfg.mPtMax = 2.0;
    cfg.mNEtaBins = 20;

    cfg.mSpecies = {211, -211, 111, 2212};

    cfg.mMaxCompositionTemplatesPerPair = 10;
    cfg.mCompositionReservoirSeed = 42;

    const std::vector<ParticleDefinition> definitions = {{211, "pi+", 0.13957039, 3}, {-211, "pi-", 0.13957039, -3}, {111, "pi0", 0.1349768, 0}, {2212, "p+", 0.93827209, 3}};

    TuneAccumulator accumulator(cfg, definitions);

    require(accumulator.processedEvents() == 0, "new accumulator should contain zero events");

    // Event 1:
    //
    // Nch = 3:
    //   pi+
    //   pi-
    //   mu-
    //
    // Nselected = 4:
    //   pi+
    //   pi-
    //   pi0
    //   proton
    //
    // The muon is intentionally not part of cfg.mSpecies.
    InputEvent event1 = {{211, 0.5, 0.1, 3}, {-211, 0.7, -0.2, -3}, {2212, 1.2, 1.0, 3}, {111, 0.4, 0.3, 0}, {13, 0.8, 0.1, -3}};

    accumulator.processEvent(event1);

    // Event 2:
    //
    // pi+ contributes to activity and composition.
    // pT=2.5 is deliberately above ptMax.
    //
    // pi0 is outside particle acceptance.
    InputEvent event2 = {{211, 2.5, 0.1, 3},
                         {111, 0.6, 2.1, 0}};

    accumulator.processEvent(event2);

    // Empty event.
    accumulator.processEvent({});

    require(accumulator.processedEvents() == 3, "expected three processed events");
    require(accumulator.activityOverflowEvents() == 0, "unexpected activity overflow");
    require(accumulator.ptOverflowParticles() == 1, "expected one pT overflow");
    require(accumulator.speciesMultiplicityOverflowEvents() == 0, "unexpected species multiplicity overflow");
    require(!accumulator.finalized(), "accumulator finalized too early");

    accumulator.finalize();

    require(accumulator.finalized(), "accumulator should be finalized");

    const auto& tune = accumulator.tune();

    require(tune.mFinalized, "Tune::finalized should be true");
    require(tune.mNEvents == 3, "expected Tune::nEvents == 3");
    require(tune.mActivityOverflowEvents == 0, "wrong activity overflow count");
    require(tune.mPtOverflowParticles == 1, "wrong pT overflow count");

    const auto* pionPlus = tune.findSpecies(211);
    const auto* pionMinus = tune.findSpecies(-211);
    const auto* pionZero = tune.findSpecies(111);
    const auto* proton = tune.findSpecies(2212);

    require(pionPlus != nullptr, "pi+ missing");
    require(pionMinus != nullptr, "pi- missing");
    require(pionZero != nullptr, "pi0 missing");
    require(proton != nullptr, "proton missing");

    require(tune.findSpecies(13) == nullptr, "muon should not be a configured species");
    requireClose(pionPlus->mHPtVsActivity.GetEntries(), 2.0, 1e-12, "wrong number of pi+ entries");
    requireClose(pionMinus->mHPtVsActivity.GetEntries(), 1.0, 1e-12, "wrong number of pi- entries");
    requireClose(proton->mHPtVsActivity.GetEntries(), 1.0, 1e-12, "wrong number of proton entries");
    requireClose(pionZero->mHPtVsActivity.GetEntries(), 1.0, 1e-12, "wrong number of pi0 entries");
    requireClose(tune.mPNch.Integral(1, tune.mPNch.GetNbinsX()), 1.0, 1e-12, "pNch is not normalized");
    requireClose(tune.mPNSelected.Integral(1, tune.mPNSelected.GetNbinsX()), 1.0, 1e-12, "pNSelected is not normalized");
    require(tune.numberOfCompositionPairs() == 3, "expected three composition pairs");
    require(tune.numberOfCompositionTemplates() == 3, "expected three composition templates");
    require(tune.compositionPairIndex(3, 4) >= 0, "missing composition pair (3,4)");
    require(tune.compositionPairIndex(1, 1) >= 0, "missing composition pair (1,1)");
    require(tune.compositionPairIndex(0, 0) >= 0, "missing composition pair (0,0)");
    // State-machine protection.
    bool threw = false;

    try {
      accumulator.processEvent(event1);
    } catch (const std::runtime_error&) {
      threw = true;
    }

    require(threw, "processEvent after finalize should throw");

    std::cout << "TuneAccumulator test passed\n";
    return 0;

  } catch (const std::exception& error) {
    std::cerr << "TuneAccumulator test failed: "
              << error.what()
              << '\n';

    return 1;
  }
}
