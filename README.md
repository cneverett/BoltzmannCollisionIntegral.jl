<img align="left" width="75" height="75" src="https://github.com/cneverett/BoltzmannCollisionIntegral.jl/blob/main/docs/src/assets/logo.svg" alt="BoltzmannCollisionIntegral.jl icon">

# BoltzmannCollisionIntegral.jl

[![][docs-latest-img]][docs-latest-url]
[![Build Status][gha-img]][gha-url]

[docs-latest-img]: https://img.shields.io/badge/Docs-Stable-lightgrey.svg
[docs-latest-url]: https://cneverett.github.io/BoltzmannCollisionIntegral.jl/dev/

[gha-img]: https://github.com/cneverett/BoltzmannCollisionIntegral.jl/actions/workflows/CI.yml/badge.svg?branch=main
[gha-url]: https://github.com/cneverett/BoltzmannCollisionIntegral.jl/actions/workflows/CI.yml?query=branch%3Amain

## Overview

`BoltzmannCollisionIntegral.jl` is a [Julia](http://julialang.org/) package for the evaluation of the relativistic Boltzmann collision integral for binary interactions $(12\rightleftharpoons34)$:
```math
    C(\boldsymbol{p}_3)=\int\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}\left[\frac{f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)}{1+\delta_{12}}W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)- \frac{f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)}{1+\delta_{12}}W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu)\right],
```
via momentum discretisation and Monte-Carlo sampling. The distribution functions $f(\boldsymbol{p})$ for the particles involved are assumed to be anisotropic (only axisymmetric is currently implemented). The collision integral is split into two components, an emission spectra dictating the rate of gain of particles emerging from a given binary interaction and an absorption spectra dictating the rate of loss of particles from the same interaction.

Evaluation is performed by assuming that the distribution function is constant over some set of discrete domains in phase space and then integrating over those domains using a Monte-Carlo method. For more information see the [Documentation](https://cneverett.github.io/BoltzmannCollisionIntegral.jl/).

Evaluation currently supports both single and multithreaded operation (with multi-CPU acceleration planned). Data is exported in the [JLD2](https://github.com/JuliaIO/JLD2.jl) file format.

## Usage
`BoltzmannCollisionIntegral.jl` is available to download from the [Julia package
manager](https://pkgdocs.julialang.org/v1/). Inside a Julia session, enter the package manager with `]`, then run the command 
Install the package using 
```julia 
pkg> add BoltzmannCollisionIntegral
```
then load the package by running
```julia
using BoltzmannCollisionIntegral
```

To perform an evaluation of the emission and absorption spectra, an example script `Run_Integration.jl` for selecting the binary interaction, discrete phase space bounds and integration conditions is located under the `src/Common/` folder of the package. It is recommended to copy this script and place it in your working folder and edit the fields as required. Then simply run
```julia-repl
include("Run_Integration.jl")
``` 

See [Getting Started](https://cneverett.github.io/BoltzmannCollisionIntegral.jl/dev/quickstart/) for in depth detail.