#=

This defines the differential/total cross section function and its normalisation for specific interactions, it therefore requires the particle names to be defined as global constants for use in determing the correct cross sections functions to compile.

=#

# Dependancies  
#include("Constants.jl")

# ====================== Hard Sphere Collisions ============================== #

"""
    dsigmadt_SphSphSphSph(sSmol,sBig,tSmol,tBig,uSmol,uBig)

returns the differential cross section for the binary interaction of hard spheres with normalised masses ``m_1,m_2,m_3,m_4=m_{\\text{Sph}}``. Normalised by ``πR_{Sph}^2``.

```math
\\frac{dσ}{dt} = \\frac{1}{s-4m_{\\text{Sph}}^2}
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2=4m_{\\text{Sph}}^2``
- `tSmol::Float32` : ``t - tBig``
- `tBig::Float32` : ``(m_3-m_1)^2=0``
- `uSmol::Float32` : ``u - uBig``
- `uBig::Float32` : ``(m_2-m_3)^2=0``
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

returns the total cross section for the binary interaction of hard spheres with normalised masses (wrt electron mass) ``m_1,m_2,m_3,m_4=m_\\text{Sph}``. Normalised by ``πR_{Sph}^2``.

```math
σ = \\frac{1}{2}
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

returns the differential cross section for electron positron annihilation to two photons. Berestetskii 1982 (88.4). Masses and momenta are normalised by the rest mass of the electron ``m_{\\text{Ele}}`` and the cross section is normalised by ``3σ_T``.

```math
\\frac{dσ}{dt} = -\\frac{1}{s(s-4)}\\left(\\left(\\frac{1}{t-1}+\\frac{1}{u-1}\\right)^2+\\left(\\frac{1}{t-1}+\\frac{1}{u-1}\\right)-\\frac{1}{4}\\left(\\frac{t-1}{u-1}+\\frac{u-1}{t-1}\\right)\\right)
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2 = 4 ∴ s = sSmol + 4``
- `tSmol::Float32` : ``t - tBig``
- `tBig::Float32` : ``(m_3-m_1)^2 = 1 ∴ t = tSmol + 1``
- `uSmol::Float32` : ``u - uBig``
- `uBig::Float32` : ``(m2-m3)^2 = 1 ∴ u = uSmol + 1``
"""
function dsigmadt_ElePosPhoPho(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)

    # -(1/(s(s-4)))*((1/(t-1)+1/(1-s-t))^2+(1/(t-1)+1/(1-s-t))-(1/4)*((t-1)/(1-s-t)+(1-s-t)/(t-1)))
    
    s = sSmol+sBig
    -(1/((s)*(sSmol)))*((1/(tSmol)+1/(uSmol))^2+(1/(tSmol)+1/(uSmol))-(1/4)*((tSmol)/(uSmol)+(uSmol)/(tSmol)))

end

"""
    sigma_ElePosPhoPho(sSmol,sBig)

returns the total cross section for electron positron annihilation to two photons. Berestetskii 1982 (88.6). Masses and momenta are normalised by the rest mass of the electron ``m_{\\text{Ele}}`` and the cross section is normalised by ``3σ_T``.

```math
σ = \\frac{1}{4s^2(s-4)}\\left((s^2+4s-8)\\log\\left(\\frac{\\sqrt{s}+\\sqrt{s-4}}{\\sqrt{s}-\\sqrt{s-4}}\\right)-(s+4)\\sqrt{s(s-4)}\\right)
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2 = 4 ∴ s = sSmol + 4``
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

returns the differential cross section for photon-photon annihilation to electron-positron pair. (Inverse proceess of electron positron annihilation to two photons). Masses and momenta are normalised by the rest mass of the electron ``m_{\\text{Ele}}`` and the cross section is normalised by ``3σ_T``.

```math
\\frac{dσ}{dt} = -\\frac{1}{s^2}\\left(\\left(\\frac{1}{t-1}+\\frac{1}{u-1}\\right)^2+\\left(\\frac{1}{t-1}+\\frac{1}{u-1}\\right)-\\frac{1}{4}\\left(\\frac{t-1}{u-1}+\\frac{u-1}{t-1}\\right)\\right)
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2 = 0 ∴ s = sSmol``
- `tSmol::Float32` : ``t - tBig``
- `tBig::Float32` : ``(m_3-m_1)^2 = 1 ∴ t = tSmol + 1``
- `uSmol::Float32` : ``u - uBig``
- `uBig::Float32` : ``(m_2-m_3)^2 = 1 ∴ u = uSmol + 1``
"""
function dsigmadt_PhoPhoElePos(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)

    # -(1/(s^2))*((1/(t-1)+1/(1-s-t))^2+(1/(t-1)+1/(1-s-t))-(1/4)*((t-1)/(1-s-t)+(1-s-t)/(t-1)))
    
    -(1/(sSmol^2))*((1/(tSmol)+1/(uSmol))^2+(1/(tSmol)+1/(uSmol))-(1/4)*((tSmol)/(uSmol)+(uSmol)/(tSmol)))

end

"""
    sigma_PhoPhoElePos(sSmol,sBig)

returns the total cross section for photon-photon annihilation to electron-positron pair. Masses and momenta are normalised by the rest mass of the electron ``m_{\\text{Ele}}`` and the cross section is normalised by ``3σ_T``.

```math
σ = \\frac{1}{2s^3}\\left((s^2+4s-8)\\log\\left(\\frac{\\sqrt(s)+\\sqrt(s-4)}{\\sqrt(s)-\\sqrt(s-4)}\\right)-(s+4)\\sqrt{s(s-4)}\\right)
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2 = 0 ∴ s = sSmol``
"""
function sigma_PhoPhoElePos(sSmol::Float32,sBig::Float32)

    #(1/(2*s^3))*((s^2+4*s-8)*log((sqrt(s)+sqrt(s-4))/(sqrt(s)-sqrt(s-4)))-(s+4)*sqrt(s*(s-4)))
    s = sSmol+sBig
    (1/(2*s^3))*((s^2+4*s-8)*log((2*s-4+2*sqrt(s*(s-4)))/(4))-(s+4)*sqrt(s*(s-4)))

end

const dsigmadtNorm_PhoPhoElePos = 3*σT;
const sigmaNorm_PhoPhoElePos = 3*σT;

# ==================================================================== # 


# =============== Electron Compton Scattering ======================== #

"""
    dsigmadt_ElePhoElePho(sSmol,sBig,tSmol,tBig,uSmol,uBig)

returns the differential cross section for electron-photon scattering (Compton) scattering. Berestetskii 1982 (86.6). Masses and momenta are normalised by the rest mass of the electron ``m_{\\text{Ele}}`` and the cross section is normalised by ``3σ_T``.

```math
\\frac{d\\sigma_{e\\gamma\\rightarrow e\\gamma}}{dt}(s,t)=\\frac{3\\sigma_Tm_e^2}{(s-m_e^2)^2}\\left[\\left(\\frac{m_e^2}{s-m_e^2}+\\frac{m_e^2}{u-m_e^2}\\right)^2+\\left(\\frac{m_e^2}{s-m_e^2}+\\frac{m_e^2}{u-m_e^2}\\right)-\\frac{1}{4}\\left(\\frac{s-m_e^2}{u-m_e^2}+\\frac{u-m_e^2}{s-m_e^2}\\right)\\right]
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2 = 1 ∴ s = sSmol + 1``
- `tSmol::Float32` : ``t - tBig``
- `tBig::Float32` : ``(m_3-m_1)^2 = 0 ∴ t = tSmol``
- `uSmol::Float32` : ``u - uBig``
- `uBig::Float32` : ``(m_2-m_3)^2 = 1 ∴ u = uSmol + 1``
"""
function dsigmadt_ElePhoElePho(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)

    # -(1/(s-1)^2)*((1/(s-1)+1/(u-1))^2+(1/(s-1)+1/(u-1))-(1/4)*((s-1)/(u-1)+(u-1)/(s-1)))
    
    (1/(sSmol)^2)*((1/(sSmol)+1/(uSmol))^2+(1/(sSmol)+1/(uSmol))-(1/4)*((sSmol)/(uSmol)+(uSmol)/(sSmol)))

end

"""
    sigma_ElePhoElePho(sSmol,sBig)

returns the total cross section for electron-photon (Compton) scattering. Berestetskii 1982 (86.16). Masses and momenta are normalised by the rest mass of the electron ``m_{\\text{Ele}}`` and the cross section is normalised by ``3σ_T``.

```math
\\sigma_{e\\gamma\\rightarrow e\\gamma}(s)=\\frac{3\\sigma_Tm_e^2}{4(s-m_e^2)}\\left[(1-\\frac{4m_e^2}{\\left(s-m_e^2\\right)}-\\frac{8m_e^4}{\\left(s-m_e^2\\right)^2})\\log\\left(1+\\frac{s-m_e^2}{m_e^2}\\right)+\\frac{1}{2}+\\frac{8m_e^2}{s-m_e^2}-\\frac{m_e^4}{2s^2}\\right]
```

# Arguments
- `sSmol::Float32` : ``s - sBig``
- `sBig::Float32` : ``(m_1+m_2)^2 = 1 ∴ s = sSmol + 1``
"""
function sigma_ElePhoElePho(sSmol::Float32,sBig::Float32)

    #(1/(4*(s-1)))*((1-4/(s-1)-8/(s-1)^2)*log(1+1/(s-1))+1/2+8/(s-1)-1/(2*s^2))
    s = sBig+sSmol
    (1/(4*(sSmol)))*((1-4/(sSmol)-8/(sSmol)^2)*log(s)+1/2+8/(sSmol)-1/(2*s^2))

end

const dsigmadtNorm_ElePhoElePho = 3*σT;
const sigmaNorm_ElePhoElePho = 3*σT;

# ==================================================================== # 

# ==================================================================== #
# ==================================================================== #



#=sSmol = 2000f0
sBig = 1f0
s = sSmol+sBig

1/s-1
2-s-1
uSmol = -1000f0


uBig = 1f0
u = uSmol+uBig
tBig = 0f0
tSmol = 2-s-u-tBig
t = tSmol+tBig

sigma_ElePhoElePho(sSmol,sBig)
dsigmadt_ElePhoElePho(sSmol,sBig,tSmol,tBig,uSmol,uBig)=#

