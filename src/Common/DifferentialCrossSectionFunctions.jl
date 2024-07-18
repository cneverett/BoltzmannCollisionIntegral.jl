#=

This defines the differential/total cross section function and its normalisation for specific interactions, it therefore requires the particle names to be defined as global constants for use in determing the correct cross sections functions to compile.

=#

# Dependancies  
include("MyPhysicalConstants.jl")

# ====================== Hard Sphere Collisions ============================== #

"""
    dsigmadt_SphSphSphSph(sSmol,sBig,tSmol,tBig,uSmol,uBig)

returns the differential cross section for the binary interaction of hard spheres with normalised masses m1,m2,m3,m4

```math
\frac{dσ}{dt} = \frac{1}{s-4μ_{\text{Sph}}^2}
```

# Arguments
- `sSmol::Float32` : s - sBig
- `sBig::Float32` : (m1+m2)^2
- `tSmol::Float32` : t - tBig
- `tBig::Float32` : (m3-m1)^2
- `uSmol::Float32` : u - uBig
- `uBig::Float32` : (m2-m3)^2
"""
function dsigmadt_SphSphSphSph(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)

    #=
        1f0/(s-4*muSph^2)
        sBig = (m3+m4)^2=4musph^2
        sSmol = s - sBig
    =#
    1f0/(sSmol) 

end

"""
    sigma_SphSphSphSph(sSmol,sBig)

returns the total cross section for the binary interaction of hard spheres with normalised masses m1,m2,m3,m4

```math
σ = \frac{1}{2}
```

# Arguments
- `sSmol::Float32` : s - sBig
- `sBig::Float32` : (m1+m2)^2
"""
function sigma_SphSphSphSph(sSmol::Float32,sBig::Float32)
    
    1f0/2f0 # factor of 2 accounts for identical final states

end

const dsigmadtNorm_SphSphSphSph = Float32(pi)*(2f0*RSph)^2
const sigmanNorm_SphSphSphSph = Float32(pi)*(2f0*RSph)^2

# ======================================================================= #

#============== Electron Positron Annihilation to Two Photons ==========#

"""
    dsigmadt_ElePosPhoPho(sSmol,sBig,tSmol,tBig,uSmol,uBig)

returns the differential cross section for electron positron annihilation to two photons. Berestetskii 1982 (88.4)

```math
\frac{dσ}{dt} = -\frac{1}{s(s-4)}\left(\left(\frac{1}{t-1}+\frac{1}{u-1}\right)^2+\left(\frac{1}{t-1}+\frac{1}{u-1}\right)-\frac{1}{4}\left(\frac{t-1}{u-1}+\frac{u-1}{t-1}\right)\right)
```

# Arguments
- `sSmol::Float32` : s - sBig
- `sBig::Float32` : (m1+m2)^2 = 4 (normalised units) -> s = sSmol + 4
- `tSmol::Float32` : t - tBig
- `tBig::Float32` : (m3-m1)^2 = 1 (normalised units) -> t = tSmol + 1
- `uSmol::Float32` : u - uBig
- `uBig::Float32` : (m2-m3)^2 = 1 (normalised units) -> u = uSmol + 1
"""
function dsigmadt_ElePosPhoPho(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)

    # -(1/(s(s-4)))*((1/(t-1)+1/(1-s-t))^2+(1/(t-1)+1/(1-s-t))-(1/4)*((t-1)/(1-s-t)+(1-s-t)/(t-1)))
    
    s = sSmol+sBig
    -(1/((s)*(sSmol)))*((1/(tSmol)+1/(uSmol))^2+(1/(tSmol)+1/(uSmol))-(1/4)*((tSmol)/(uSmol)+(uSmol)/(tSmol)))

end

"""
    sigma_ElePosPhoPho(sSmol,sBig)

returns the total cross section for electron positron annihilation to two photons. Berestetskii 1982 (88.6)

```math
σ = \frac{1}{4s^2(s-4)}\left((s^2+4s-8)\log\left(\frac{\sqrt{s}+\sqrt{s-4}}{\sqrt{s}-\sqrt{s-4}}\right)-(s+4)\sqrt{s(s-4)}\right)
```

# Arguments
- `sSmol::Float32` : s - sBig
- `sBig::Float32` : (m1+m2)^2 = 4 (normalised units) -> s = sSmol + 4
"""
function sigma_ElePosPhoPho(sSmol::Float32,sBig::Float32)

    #(1/(4*s^2*(s-4)))*((s^2+4*s-8)*log((sqrt(s)+sqrt(s-4))/(sqrt(s)-sqrt(s-4)))-(s+4)*sqrt(s*(s-4)))

    s = sSmol+sBig
    (1/(4*(sSmol)*s^2))*((sSmol^2+12*sSmol+24)*log((s+sSmol+2*sqrt(sSmol*s))/(sBig))-(sSmol+8)*sqrt((s)*(sSmol)))

end

const dsigmadtNorm_ElePosPhoPho = 3*σT;
const sigmaNorm_ElePosPhoPho = 3*σT;

# ==================================================================== # 


#======== Electron Positron Pair Production from Two Photons ==========#

"""
    dsigmadt_PhoPhoElePos(sSmol,sBig,tSmol,tBig,uSmol,uBig)

returns the differential cross section for photon-photon annihilation to electron-positron pair. (Inverse proceess of electron positron annihilation to two photons) 

```maths
\frac{dσ}{dt} = -\frac{1}{s^2}\left(\left(\frac{1}{t-1}+\frac{1}{u-1}\right)^2+\left(\frac{1}{t-1}+\frac{1}{u-1}\right)-\frac{1}{4}\left(\frac{t-1}{u-1}+\frac{u-1}{t-1}\right)\right)
```

# Arguments
- `sSmol::Float32` : s - sBig
- `sBig::Float32` : (m1+m2)^2 = 0 (normalised units) -> s = sSmol
- `tSmol::Float32` : t - tBig
- `tBig::Float32` : (m3-m1)^2 = 1 (normalised units) -> t = tSmol + 1
- `uSmol::Float32` : u - uBig
- `uBig::Float32` : (m2-m3)^2 = 1 (normalised units) -> u = uSmol + 1
"""
function dsigmadt_PhoPhoElePos(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)

    # -(1/(s^2))*((1/(t-1)+1/(1-s-t))^2+(1/(t-1)+1/(1-s-t))-(1/4)*((t-1)/(1-s-t)+(1-s-t)/(t-1)))
    
    -(1/(sSmol^2))*((1/(tSmol)+1/(uSmol))^2+(1/(tSmol)+1/(uSmol))-(1/4)*((tSmol)/(uSmol)+(uSmol)/(tSmol)))

end

"""
    sigma_PhoPhoElePos(sSmol,sBig)

returns the total cross section for photon-photon annihilation to electron-positron pair.

```math
sigma = \frac{1}{2s^3}\left((s^2+4s-8)\log\left(\frac{\sqrt(s)+\sqrt(s-4)}{\sqrt(s)-\sqrt(s-4)}\right)-(s+4)\sqrt{s(s-4)}\right)
```

# Arguments
- `sSmol::Float32` : s - sBig
- `sBig::Float32` : (m1+m2)^2 = 0 (normalised units) -> s = sSmol
"""
function sigma_PhoPhoElePos(sSmol::Float32,sBig::Float32)

    #(1/(2*s^3))*((s^2+4*s-8)*log((sqrt(s)+sqrt(s-4))/(sqrt(s)-sqrt(s-4)))-(s+4)*sqrt(s*(s-4)))
    s = sSmol+sBig
    (1/(2*s^3))*((s^2+4*s-8)*log((2*s-4+2*sqrt(s*(s-4)))/(4))-(s+4)*sqrt(s*(s-4)))

end

const dsigmadtNorm_PhoPhoElePos = 3*σT;
const sigmaNorm_PhoPhoElePos = 3*σT;

# ==================================================================== # 

# ==================================================================== #
# ==================================================================== #
