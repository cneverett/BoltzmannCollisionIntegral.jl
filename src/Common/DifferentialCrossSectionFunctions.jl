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

#=

sSmol = 2000f0
sBig = 1f0
s = sSmol+sBig

1/s-1
2-s-1
uSmol = -1000f0

1.892324e6#3.7846492e6#1.0#-3.784648e6#0.0#0.8574219#1.0#Float32[4239.24, 0.78002834, 0.37422585]#Float32[262.28375, -0.28028643, 1.5766257]#Float32[4239.385, 0.7800194, 0.37416828]

425440.9#850878.4#1.0#-850881.94#0.0#3.8998935#1.0#Float32[5719.873, 0.18377423, 1.0450137]#Float32[5719.476, 0.18395984, 1.0450138]#Float32[91.98793, -0.23729122, 0.6241733]

p3v = Float32[4239.24, 0.78002834, 0.37422585]

p4v = Float32[5719.873, 0.18377423, 1.0450137]

p1v = Float32[5719.476, 0.18395984, 1.0450138]
p2v = Float32[91.98793, -0.23729122, 0.6241733]

sSmol::Float32 = 850878.4
sBig::Float32 = 1.0
tSmol::Float32 = -850881.9
tBig::Float32 = 0.0
uSmol::Float32 = 3.8998935
uBig::Float32 = 1.0

m3 = 1f0
m2 = 0f0
m1 = 1f0
m4 = 0f0

p1 = p1v[1]
p2 = p2v[1]
ct1 = p1v[2]
ct2 = p2v[2]
st1 = sqrt(1f0-ct1^2)
st2 = sqrt(1f0-ct2^2)
Es2::Float32 = m2 != 0f0 ? (p2^2)/(sqrt(m2^2+p2^2)+m2) : p2
Es2s::Float32 = Es2/p2
Es1::Float32 = m1 != 0f0 ? (p1^2)/(sqrt(m1^2+p1^2)+m1) : p1
Es1s::Float32 = Es1/p1
E1::Float32 = m1 + Es1
E2::Float32 = m2 + Es2

p3 = p3v[1]
ct3 = p3v[2]
st3 = sqrt(1f0-ct3^2)
ch3h1 = cospi(p3v[3]-p1v[3])
ch3h2 = cospi(p3v[3]-p2v[3])
Es3::Float32 = m3 != 0f0 ? (p3^2)/(sqrt(m3^2+p3^2)+m3) : p3
Es3s::Float32 = Es3/p3
E3::Float32 = m3 + Es3
tSmol::Float32 = -2*(m1*Es3 + m3*Es1 + p3*p1*(Es3s*Es1s-(ct3*ct1+ch3h1*st3*st1))) 
uSmol::Float32 = -2*(m3*Es2 + m2*Es3 + p2*p3*(Es2s*Es3s-(ct2*ct3+ch3h2*st2*st3)))

p4::Float32 = p4v[1]
ct4::Float32 = p4v[2] # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
st4::Float32 = sqrt(1f0-p4v[2]^2)
ch4h1::Float32 = cospi(p4v[3]-p1v[3])
ch4h2::Float32 = cospi(p4v[3]-p2v[3])
Es4::Float32 = m4 != 0f0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
Es4s::Float32 = Es4/p4
E4::Float32 = Es4 + m4
uSmol::Float32 = -2*(m1*Es4 + m4*Es1 + p4*p1*(p4/(sqrt(m4^2+p4^2)+m4)*p1/(sqrt(m1^2+p1^2)+m1) - (ct4*ct1+ch4h1*st4*st1)))
tSmol::Float32 = -2*(m4*Es2 + m2*Es4 + p2*p4*(p4/(sqrt(m4^2+p4^2)+m4)*p2/(sqrt(m2^2+p2^2)+m2) - (ct2*ct4+ch4h2*st2*st4)))

uSmol::Float32 = -2*p4*p1*( - ct4*ct1 - ch4h1*st4*st1 + Es4s*Es1s + m1*Es4s/p1 + m4*Es1s/p4 )
tSmol::Float32 = -2*p2*p4*(- ct2*ct4 - ch4h2*st2*st4 + Es2s*Es4s + m4*Es2s/p4 + m2*Es4s/p2 )

p4*p1
m1*Es4s/p1 
m4*Es1s/p4 
Es4s*Es1s 
-ct4*ct1
-ch4h1*st4*st1

( - ct4*ct1 - ch4h1*st4*st1 + Es4s*Es1s + m1*Es4s/p1 + m4*Es1s/p4 )

s = sBig+sSmol
t = tBig+tSmol
u = uBig+uSmol
s+t+u

sigma_ElePhoElePho(sSmol,sBig)
dsigmadt_ElePhoElePho(sSmol,sBig,tSmol,tBig,-5.9f0,uBig)
dsigmadt_ElePhoElePho(sSmol,sBig,tSmol,tBig,uSmol,uBig)

p3v_test = (copy(p3v))
p3pv_test = copy(p3v_test)
Momentum3Value!(p3v_test,p3pv_test,p1v,p2v,m1,m2,m3,m4)
p3v_test
p3pv_test

p1v64 = copy(Float64.(p1v))
p2v64 = copy(Float64.(p2v))
p3v64 = copy(Float64.(p3v))
p3pv64 = copy(p3v64)
Momentum3Value!(p3v64,p3pv64,p1v64,p2v64,m1,m2,m3,m4)
p3v64
p3pv64

p4v64 = copy(Float64.(p4v))
SValue4(p4v64,p1v64,p2v64,dsigmadt_ElePhoElePho64,m1,m2,m3,m4)
SValue4(p4v,p1v,p2v,dsigmadt_ElePhoElePho,m1,m2,m3,m4)


function Momentum3Value!(p3v::Vector{Float64},p3pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},mu1::Float32,mu2::Float32,mu3::Float32,mu4::Float32)

    # set normalised masses (defined in Init.jl)
    m1::Float64 = mu1
    m2::Float64 = mu2
    m3::Float64 = mu3
    m4::Float64 = mu4 

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    ct3::Float64 = p3v[2] 
    ct1::Float64 = p1v[2]
    ct2::Float64 = p2v[2] 

    st3::Float64 = sqrt(1f0-ct3^2)
    st1::Float64 = sqrt(1f0-ct1^2)
    st2::Float64 = sqrt(1f0-ct2^2)

    ch1h3::Float64 = cospi(p3v[3]-p1v[3])
    ch1h4::Float64 = cospi(p3v[3]-p2v[3])
    ch3h4::Float64 = cospi(p1v[3]-p2v[3])

    m32::Float64 = m3^2
    m42::Float64 = m4^2
    m12::Float64 = m1^2
    m22::Float64 = m2^2

    p12::Float64 = p1^2
    p22::Float64 = p2^2

    sqm1p1::Float64 = sqrt(m12+p12)
    sqm2p2::Float64 = sqrt(m22+p22)

    ct3ct1::Float64 = ct3*ct1 
    ct3ct2::Float64 = ct3*ct2
    ct1ct2::Float64 = ct1*ct2
    st3st1::Float64 = st3*st1
    st3st2::Float64 = st3*st2
    st1st2::Float64 = st1*st2

    A1::Float64 = m1 != 0f0 ? p12/(sqm1p1+m1) : p1
    A2::Float64 = m2 != 0f0 ? p22/(sqm2p2+m2) : p2

    #sqm1p1sqm2p2 = sqm1p1*sqm2p2

    #p1p2 = p1*p2

    # reset p3v values
    p3v[1] = 0f0 
    p3pv[1] = 0f0

    C3sqr::Float64 = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))^2*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*(A2-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2)) 

    C2::Float64 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch3h4*st1st2))

    C4::Float64 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))

    # C3sqr == 0 was causing issues with SValue calculation often leading to deltacorrect = 0 so we are going to ignore this point and tread it as if p3 were complex.
    #=if C3sqr == 0 # only one state and p3 cannont equal zero

        NumStates = 1
        p3_physical = false
        p3p_physical = false

        p3 = C2/C4
        if p3 == 0f0
            NumStates = 0
        else
            if p3 > 0f0
                p3v[1] = p3
            else
                p3v[1] = -p3
                p3v[2] *= -1
                p3v[3] = mod(p3v[3]+1f0,2f0)
            end
            if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
                p3_physical = true
            end
        end

    else=#if C3sqr > 0

        NumStates = 2

        p3_physical = false
        p3p_physical = false

        C3 = 4*sqrt(C3sqr)
        p3 = (C2-C3)/C4

        if p3 == 0f0
            NumStates = 1
        else
            if p3 > 0f0
                p3v[1] = p3
            else
                p3v[1] = -p3
                p3v[2] *= -1
                p3v[3] = mod(p3v[3]+1f0,2f0)
            end

            if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
                p3_physical = true
            end
        end

        if NumStates == 2
            p3p = (C2+C3)/C4
            if p3p == 0f0
                NumStates = 1
            else
                if p3p > 0f0
                    p3pv[1] = p3p
                    
                else
                    p3pv[1] = -p3p
                    p3pv[2] *= -1
                    p3pv[3] = mod(p3pv[3]+1f0,2f0)
                end

                if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3p^2/(sqrt(m32+p3p^2)+m3)) > m3-m1-m2
                    p3p_physical = true
                end
            end
        else # NumStates == 1
            p3 = (C2+C3)/C4
            if p3 == 0
                NumStates = 0
            else
                if p3 > 0f0
                    p3v[1] = p3
                else
                    p3v[1] = -p3
                    p3v[2] *= -1
                    p3v[3] = mod(p3v[3]+1f0,2f0)
                end
                if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
                    p3_physical = true
                end
            end    

        end

    else # imaginary C3sqr < 0f0

        NumStates = 1 # two states but both in same bin so same as one
        p3p_physical = false
        p3_physical = false

        p3Real = C2/C4
        if p3Real == 0f0
            NumStates = 0
        else
            if p3Real < 0f0
            p3v[2] *= -1
            p3v[3] = mod(p3v[3]+1f0,2f0)
            #p3pv[2] *= -1
            #p3pv[3] = mod(p3pv[3]+1f0,2f0)
            end
        end

    end

    return p3_physical, p3p_physical, NumStates

end

function SValue4(p4v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},dsigmadt::Function,mu1::Float32,mu2::Float32,mu3::Float32,mu4::Float32)

    # define normalise masses
    m1 = mu1
    m2 = mu2 
    m3 = mu3
    m4 = mu4 

    # pre-defining terms for efficiency 
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    ct1::Float64 = p1v[2] #cospi(p1v[2])
    ct2::Float64 = p2v[2] #cospi(p2v[2]) 

    st1::Float64 = sqrt(1f0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float64 = sqrt(1f0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2

    Es1::Float64 = m1 != 0f0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    E1::Float64 = Es1 + m1
 
    Es2::Float64 = m2 != 0f0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    E2::Float64 = Es2 + m2

    sBig::Float64 = (m1+m2)^2
    sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float64 = (1/E1)*(1/E2)*(2*InvarientFlux2Small(sSmol,m1,m2))

    p4::Float64 = p4v[1]
    ct4::Float64 = p4v[2] # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    st4::Float64 = sqrt(1f0-p4v[2]^2)
    ch4h1::Float64 = cospi(p4v[3]-p1v[3])
    ch4h2::Float64 = cospi(p4v[3]-p2v[3])
    Es4::Float64 = m4 != 0f0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
    E4::Float64 = Es4 + m4

    # u = uBig + uSmol
    uBig::Float64 = (m4-m1)^4
    uSmol::Float64 = -2*(m1*Es4 + m4*Es1 + Es4*Es1 - p4*p1*(ct4*ct1+ch4h1*st4*st1))
    # t = tBig + tSmol
    tBig::Float64 = (m2-m4)^2
    #tSmol::Float64 = m12+m22+m32+m42 - sBig - uBig - tBig - sSmol - uSmol # Leads to Float64 issues, better to calculate directly 
    tSmol::Float64 = -2*(m4*Es2 + m2*Es4 + Es2*Es4 - p2*p4*(ct2*ct4+ch4h2*st2*st4))

    deltacorrect::Float64 = (Es1*p4 - Es4*p1*(ct4*ct1+ch4h1*st4*st1) + Es2*p4 - Es4*p2*(ct4*ct2+ch4h2*st4*st2)) + (m1*p4 - m4*p1*(ct4*ct1+ch4h1*st4*st1) + m2*p4 - m4*p2*(ct4*ct2+ch4h2*st4*st2))
    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    Sval = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p4^2/(deltacorrect*sign(deltacorrect)))

    if (Sval==Inf || Sval == -Inf || Sval < 0f0)
        error("ST1 Inf#$deltacorrect#$sSmol#$sBig#$tSmol#$tBig#$uSmol#$uBig#$p4v#$p1v#$p2v")  
    end

    return Sval

end

function dsigmadt_ElePhoElePho64(sSmol::Float64,sBig::Float64,tSmol::Float64,tBig::Float64,uSmol::Float64,uBig::Float64)

    # -(1/(s-1)^2)*((1/(s-1)+1/(u-1))^2+(1/(s-1)+1/(u-1))-(1/4)*((s-1)/(u-1)+(u-1)/(s-1)))
    
    (1/(sSmol)^2)*((1/(sSmol)+1/(uSmol))^2+(1/(sSmol)+1/(uSmol))-(1/4)*((sSmol)/(uSmol)+(uSmol)/(sSmol)))

end

function InvarientFlux2Small(sSmol::Float64,mu1::Float32,mu2::Float32)
    # Better accuracy for small s

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (sSmol*(sSmol+4*mu1*mu2))/4

end

=#

