# Ditto PYTHIA tuner

This is the first teacher-based tuning stage for Ditto.

## What it does

`Ditto::Tuner` runs PYTHIA minimum-bias events, accumulates the
statistics in memory, computes normalized probability tables only after the
requested number of successful PYTHIA events has been generated, and then
writes one ROOT tune file.

The complete PYTHIA configuration is read from `cfg.pythiaCard`. For pp, a
card can contain:

```text
Beams:idA = 2212
Beams:idB = 2212
SoftQCD:inelastic = on
```

with `Beams:eCM = 13000.`.

For PbPb, put the Pb-208 nuclear IDs and the Angantyr settings in the card.
This keeps the card as the single source of PYTHIA configuration.

The tune is always written as `Ditto_<card-stem>.root`, making the card
and tune unambiguous.

## Tune contents

The tune file contains:

```text
Metadata/
  nEvents
  activityEtaMax
  particleEtaMax
  ActivityEdges
  pythiaCard
  SpeciesPDG

Event/
  hNch
  pNch
  hNSelected
  pNSelected
  hActivity
  pActivity

Species/
  pdg_211/
    hCountVsActivity
    pCountGivenActivity
    hPtVsActivity
    pPtGivenActivity
    hEtaVsActivity
    pEtaGivenActivity

  pdg_-211/
    ...

  ...
```

The event activity is currently

```text
Nch = number of final charged PYTHIA particles in |eta| < activityEtaMax.
```

For each activity class, the three probability tables are normalized so that
each regular Y slice sums to one:

```text
P(N_species | Nch class)
P(pT        | species, Nch class)
P(eta       | species, Nch class)
```

Zero-particle events for a given species are explicitly filled into
`pCountGivenActivity`, so that it really is a multiplicity probability and not
a distribution conditional on the species being present.

## Particle definition

The tuner uses `Pythia8::Particle::isFinal()`.

There is no special stable/transport flag in the tuner. The tune learns the
leaf/final event record produced by the configured PYTHIA teacher. Therefore
any decay-policy changes should be made as ordinary PYTHIA settings in
the PYTHIA card.

The default species list is taken directly from `Ditto::defaultSpecies()`,
so it remains synchronized with the generator.

## ROOT usage

```cpp
.L ../Ditto.cxx+
.L DittoTuner.cxx+
.x exampleTuner.C(1000000)
```

Run these commands from the `Tuning` directory, in this order. Processing
`exampleTuner.C` by itself only compiles the wrapper macro; it does not load
the `Ditto` or `Tuner` implementations.

For serious tuning, millions of pp events are preferable, especially for
high-multiplicity and high-pT conditional bins.

## Important limitation of v1

The tune currently factorizes kinematics as

```text
P(Nch)
P(N_species | Nch class)
P(pT | species, Nch class)
P(eta | species, Nch class)
```

It does **not** yet encode the correlation between pT and eta, nor
particle-particle Delta-eta / Delta-phi correlations.

That is intentional. The next step should be to plug these tables into
Ditto, check closure against PYTHIA, and add dimensions only where the
closure test shows they are needed.
