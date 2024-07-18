# Cross Sections

Below is a list of the currently available binary interaction cross sections that have been implemented:
- Collision of hard spheres `SphSphSphSph`
- Photon pair production from electron positron annihilation `ElePosPhoPho`
- Electron positron pair production from photon pair annihilation `PhoPhoElePos`

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