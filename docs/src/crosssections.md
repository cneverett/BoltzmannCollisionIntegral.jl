# Implimented Particles and Binary Interactions

## [Particles](@id particle_list) 
Below is a table of the currently implimetned particles (i.e. their particle properties are defined within the code)
| Particle | Abr. String | Notes  | 
|----------|-------------|--------|
| Sphere   | `"Sph"`     |  Mass taken to be the mass of the Proton   | 
| Electron | `"Ele"`     |          | 
| Positron | `"Pos"`     |          | 
| Proton   | `"Pro"`     |          |

## [Implimented Interactions](@id implemented_interactions):

These binary interactions have currently been implemented:
- Collision of hard spheres `SphSphSphSph`
- Photon pair production from electron positron annihilation `ElePosPhoPho`
- Electron positron pair production from photon pair annihilation `PhoPhoElePos`

## Differential and total cross section functions:

All cross sections are to be defined in terms of the Mandelstram variables $s=(p_1^\mu+p_2^\mu)^2$, $t=(p_1^\mu-p_3^\mu)^2$ and $u=(p_2^\mu-p^3^\mu)^2. To maintain accuracy of cross sections and avoid DivZero issues when momenta is small compared to the mass of the particles (at Float32 precision), each Mandelstram variable is split into two components e.g. $s=sSmol+sBig$ where $sBig = (m_1+m_2)^2$. The latter part typically cancels with terms in the cross sections, leading to better accuracy. 

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