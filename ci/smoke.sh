#!/usr/bin/env bash

set -euo pipefail

export PYTHIA_ROOT="${CONDA_PREFIX}"

echo "=== Environment ==="
root-config --version
pythia8-config --version || true
cmake --version
c++ --version

echo "=== Configure/build/install ==="
make install JOBS=2

echo "=== Check shared-library dependencies ==="
ldd lib/libDitto.so

if ! ldd lib/libDitto.so | grep -q libcurl; then
  echo "ERROR: libDitto.so is not linked against libcurl"
  exit 1
fi

echo "=== Tiny PYTHIA tuning run ==="

root -l -b -q \
  'Tuning/exampleTuner.C(500,"Tuning/cards/pythia8_inel_136tev.cfg")'

TUNE="Tuning/tunes/Ditto_tune_pythia8_inel_136tev.root"

if [[ ! -s "${TUNE}" ]]; then
  echo "ERROR: tune was not produced"
  exit 1
fi

echo "=== Tiny Ditto generation run ==="

root -l -b -q \
  "Generation/example.C(1000,\"${TUNE}\",false)"

echo "=== CI smoke test passed ==="
