# Evolutionary game of direct reciprocity in group-structured populations

Simulation codes for evolutionary games of direct reciprocity in memory-1 strategy space in group-structured populations.

## Build

Lapack and cmake is necessary as a prerequisite.

```shell
$ brew install lapack
$ brew install cmake
```

Clone this repository including **submodules**. (Don't forget `--recursive`.)

```shell
$ git clone --recursive https://github.com/yohm/sim_grouped_direct_reciprocity.git
```

Make a build directory, and build with cmake.

```shell
$ make build
$ cd build
$ cmake ..
$ make
```

## Executables

- main_mem1_low_mut_MC
    - MC simulation for the low mutation rate limit
- main_mem1_finite_mut
    - MC simulation for finite mutation rate
- main_mem1_finite_mut_analysis
    - Solving ordinary differential equations for finite mutation rate limit

The other binaries (`test_....`) are unit-test programs.


