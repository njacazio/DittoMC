# Ditto

Fast generator for minimum-bias events.

![Ditto logo](.logo.png)

## Requirements

Ditto requires:

* [ROOT](https://root.cern/)
* [Pythia 8](https://pythia.org/)

ROOT must be available in the environment used for configuration and execution.

If Pythia 8 is installed in a non-standard location, its installation prefix can be passed through the `PYTHIA_ROOT` environment variable:

```sh
export PYTHIA_ROOT=/path/to/pythia8
```

## Build

Ditto uses CMake internally and provides a Makefile interface for the most common operations.

The default build type is `Release`, using:

```text
-O3 -DNDEBUG
```

To configure and build:

```sh
make
```

By default, the CMake build directory is:

```text
build/
```

A different build directory can be selected with:

```sh
make BUILD_DIR=mybuild
```

### Parallel build

To build using a specific number of parallel jobs:

```sh
make JOBS=8
```

### Installation

To build and install Ditto:

```sh
make install
```

The project installs its headers and libraries locally into the corresponding `include/` and `lib/` directories.

ROOT dictionary PCM files are installed in `lib/` alongside the shared libraries.

### Production build

For optimized production use:

```sh
make production
```

This builds and installs Ditto using:

```text
-O3 -DNDEBUG -march=native -mtune=native
```

These options optimize the generated code for the architecture of the machine performing the compilation.

### Custom build configuration

The CMake build type can be changed, for example:

```sh
make BUILD_TYPE=Debug
```

The release compiler flags can also be overridden:

```sh
make CXX_FLAGS_RELEASE="-O2 -DNDEBUG"
```

## Cleaning and rebuilding

To remove the build products and locally installed headers and libraries:

```sh
make clean
```

To clean, rebuild, and reinstall the project:

```sh
make rebuild
```

## Examples and validation

Several convenience targets are provided.

Run the example event generator:

```sh
make gen
```

This builds and installs Ditto before running:

```text
Generation/example.C
```

Run the tuning example:

```sh
make tune
```

which executes:

```text
Tuning/exampleTuner.C
```

Run the standard validation:

```sh
make val
```

which executes:

```text
Validation/ValidateDitto.C
```

Run the tuning validation:

```sh
make valTune
```

which executes:

```text
ValidationTune/ValidateDittoTune.C
```

## Documentation

API and developer documentation can be generated with Doxygen:

```sh
make doc