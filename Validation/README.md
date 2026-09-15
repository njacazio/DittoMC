# Ditto validation against ALICE 13 TeV identified spectra

`ValidateDitto.C` reads the Ditto output tree

```
T
└── Particles : TClonesArray<TParticle>
```

and compares generated charged pions, kaons and protons with ALICE inelastic pp data at sqrt(s)=13 TeV, |y|<0.5, from HEPData record `ins1797443` / DOI `10.17182/hepdata.100303`:

- Table 1: pi+ + pi-
- Table 2: K+ + K-
- Table 6: p + pbar

The MC histograms use the exact HEPData pT binning and are normalized as

```
(1/N_INEL) d^2N/(dpT dy)
```

with Delta y = 1 for |y|<0.5.

## Run

```cpp
root -l
.L ValidateDitto.C+
ValidateDitto("Ditto.root");
```

On first use, the macro downloads the three HEPData ROOT tables with `curl` into `./hepdata/`.

If they are already downloaded:

```cpp
ValidateDitto("Ditto.root", "hepdata",
                  "Ditto_validation.root",
                  "Ditto_vs_ALICE_13TeV.pdf",
                  false);
```

Outputs:

- `Ditto_validation.root`: MC spectra, HEPData objects and MC/data ratios.
- `Ditto_vs_ALICE_13TeV.pdf`: one comparison page for each species.

The HEPData tables exclude the common 2.6% normalization uncertainty; the macro does not add it back.
