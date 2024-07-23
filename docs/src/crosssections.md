# Implemented Particles and Binary Interactions

## Particles
Below is a table of the currently implemented particles (i.e. their particle properties are defined within the code)

| Particle | Abr. String | Notes                                        | 
| -------- | ----------- | -------------------------------------------- |
| Sphere   | `"Sph"`     |  Mass taken to be the mass of the Proton     | 
| Electron | `"Ele"`     |                                              | 
| Positron | `"Pos"`     |                                              | 
| Proton   | `"Pro"`     |                                              |

## Implemented Interactions

These binary interactions have currently been implemented:
- Collision of hard spheres `SphSphSphSph`
    - functions: [`BoltzmannCollisionIntegral.dsigmadt_SphSphSphSph`](@ref) [`BoltzmannCollisionIntegral.sigma_SphSphSphSph`](@ref)
- Photon pair production from electron positron annihilation `ElePosPhoPho`
    - functions: [`BoltzmannCollisionIntegral.dsigmadt_ElePosPhoPho`](@ref) [`BoltzmannCollisionIntegral.sigma_ElePosPhoPho`](@ref)
- Electron positron pair production from photon pair annihilation `PhoPhoElePos`
    - functions: [`BoltzmannCollisionIntegral.dsigmadt_PhoPhoElePos`](@ref) [`BoltzmannCollisionIntegral.sigma_PhoPhoElePos`](@ref)

## Adding User Defined Interactions

Users may add their own binary interaction cross sections to the `/src/commom/DifferentialCrossSectionFunctions.jl` file. Functions should be named in the following format: `"sigma_name1name2name3name4"` and `"dsigmadt_name1name2name3name4"` where the names are three letter abbreviations of the particles involved (see [Particles](#particles) for examples). The named pairs `name1name2` and `name3name4` should be in alphabetical order. Both the total cross section and differential cross sections must be provided for a single interaction.    

All cross sections are to be defined in terms of the Mandelstram variables $s=(p_1^\mu+p_2^\mu)^2$, $t=(p_1^\mu-p_3^\mu)^2$ and $u=(p_2^\mu-p_3^\mu)^2$. To maintain accuracy of cross sections and avoid DivZero issues when momenta is small compared to the mass of the particles (at `Float32` precision), each Mandelstram variable in the cross sections should be split into two components:
- `s=sSmol+sBig` where `sBig = (m_1+m_2)^2`
- `t=tSmol+tBig` where `tBig = (m_1-m_3)^2`
- `u=uSmol+uBig` where `uBig = (m_2-m_3)^2`
The "Big" part typically cancels with terms in the cross sections, leading to better accuracy. Therefore the function should defined as follows:
```julia
function sigma_name1name2name3name4(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
    ... 
end

function dsigmadt_name1name2name3name4(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
    ... 
end
```
where all of `sSmol`, `sBig`, `tSmol`, `tBig`, `uSmol`, and `uBig` must be included in the function definition irrespective of if they actually appear in the cross section. It is also important to define a normalisation for the cross sections to avoid emission and absorption terms becoming small compared to the `Float32` minimum.

## Differential and Total Cross Section Functions

!!! warning
    To ensure greater computational accuracy and prevent underflow of ``Float32`` precision values, all cross sections have a normalisation defined in the function documentation.

```@meta
CurrentModule = BoltzmannCollisionIntegral
using BoltzmannCollisionIntegral
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