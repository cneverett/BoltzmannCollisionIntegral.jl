# BinaryInteractionSpectra.jl

[![Build Status](https://github.com/cneverett/BinaryInteractionSpectra.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/cneverett/BinaryInteractionSpectra.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Overview

`BinaryInteractionSpectra.jl` is a [Julia](http://julialang.org/) package for the evaluation of the relativistic collision integral for binary interactions ``(12\\rightleftharpoons 34)``:
```math
    C(\boldsymbol{p}_1)=\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}\left[\frac{f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)}{1+\delta_{34}}W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu)- \frac{f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)}{1+\delta_{34}}W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)\right],
```
when the distribution function for the particles involved may be either isotropic (not implemented) or axisymmetric. The collision integral is split into two components, an emission spectra dictating the rate of gain of particles emerging from a given binary interaction and an absorption spectra dictating the rate of loss of particles from the same interaction.

Evaluation is performed by assuming that the distribution function is constant over some set of discrete domains in phase space and then integrating over those domains using a Monte-Carlo method. For more information see the [Documentation](https://cneverett.github.io/BinaryInteractionSpectra.jl/).

Evaluation currently supports both single and multithreaded operation.

## Usage

After installing the package, run 
```julia
using BinaryInteractionSpectra
```

To perform an evaluation of the emission and absorption spectra, an example script `Run_BinaryInteractionSpectra.jl` for selecting the binary interaction, discrete phase space bounds and integration conditions is located under the `src/Common/` folder of the package. It is recommended to copy this script and place it in your working folder and edit the fields as required. Then simply run
```julia-repl
include("Run_BinaryInteractionSpectra.jl")
``` 

See [Getting Started](https://cneverett.github.io/BinaryInteractionSpectra.jl/dev/quickstart/) for in depth detail.