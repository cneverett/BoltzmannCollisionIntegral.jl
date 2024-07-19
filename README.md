# BinaryInteractionSpectra.jl

[![Build Status](https://github.com/cneverett/BinaryInteractionSpectra.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/cneverett/BinaryInteractionSpectra.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Overview
---

`BinaryInteractionSpectra.jl` is a [Julia](http://julialang.org/) package for the evaluation of the relativistic collision integral for binary interactions (12->34) when the distribution function for the particles involved may be either isotropic (not implemented) or axisymmetric. The collision integral is split into two components, an emission spectra dictating the rate of gain of particles emerging from a given binary interaction and an absorption spectra dictating the rate of loss of particles from the same interaction.

Evaluation is performed by assuming that the distribution function is constant over some set of discrete domains in phase space and then integrating over those domains using a Monte-Carlo method. See documentation (to be added) for more detailed description.

Evaluation currently supports both single and multithreaded operation.

For more information see the [Documentation](https://cneverett.github.io/BinaryInteractionSpectra.jl/).

## Usage
---

See [Getting Started](https://cneverett.github.io/BinaryInteractionSpectra.jl/dev/quickstart/) for in depth detail. But in brief:

After installing the package, run 
```julia
using BinaryInteractionSpectra
```

To perform an evaluation of the emission and absorption spectra, an example script `Run_BinaryInteractionSpectra.jl` for selecting the binary interaction, discrete phase space bounds and integration conditions is located under the `src/Common/` folder of the package. It is recommended to copy this script and place it in your working folder and edit the fields as you require. Then simply run
```julia-repl
include("Run_BinaryInteractionSpectra.jl")
``` 