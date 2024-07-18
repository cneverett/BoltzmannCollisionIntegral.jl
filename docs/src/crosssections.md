# Available Binary Interactions

Below is a list of the currently available binary interactions:
- Collision of hard spheres
- Photon pair production from electron positron annihilation
- Electron positron pair production from photon pair annihilation

## Differential and total cross section functions:

```@meta
CurrentModule = BinaryInteractionSpectra
using BinaryInteractionSpectra
end
```

```@docs
sigma_SphSphSphSph
dsigmadt_SphSphSphSph

sigma_ElePosPhoPho
dsigmadt_ElePosPhoPho

sigma_PhoPhoElePos
dsigmadt_PhoPhoElePos
```