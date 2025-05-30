#=
Contains functions for generating random points
=#

"""
    RPointLogMomentum!(pu,pl,pv,num)

Edits the first element of `pv` with a random real-space momentum value between ``10^{pl}`` and ``10^{pu}``. This sample is chosen by first randomly picking a momentum bin in the range `1:num` and then uniformly sampling a momentum point in real-space (rather than log10 space) between pl and pu which are the momentum values at start and end of that bin. Sampling is done such there will be a constant number of points per momentum-space volume. As the momentum space between ``10^{pl}`` and ``10^{pu}`` it is a spherical shell hence the correct sampling is ``p = (U*(10^{pu})^3+(1-U)*(10^{pl})^3)^{1/3}`` with uniform ``U ∈ [0~~1]``.

Assumes ``f(x,p,μ)=f(x,\\vec{p})*(2πp^2)=const`` in bin, therefore momentum space volume element is ``\\mathrm{d}p`` and as such uniform sampling corresponds to ``U*10^{u}+(1-U)*10^{l}`` where ``U`` is a uniform random number between 0 and 1.

If instead ``f(x,\\vec{p})=const`` in bin, momentum space volume element is ``p^2 \\mathrm{d}p`` and uniform sampling corresponds to ``(10^pu)*\\sqrt[3]{U+(1-U)*10^{3pl-3pu}}`` where ``U`` is a uniform random number between 0 and 1.
"""
function RPointLogMomentum!(pv::Vector{Float64},pu::Float64,pl::Float64,num::Int64) 
    # Inputs a momentum vector and momentum bounds and mutates first of said vector
    bin = rand(1:num)
    l = (pl + (pu-pl)*(bin-1)/num)
    u = (pl + (pu-pl)*(bin)/num)

    U = rand(Float64)
    #pv[1] = (10^u)*cbrt(U+(1-U)*1f3^(l-u)) 
    pv[1] = U*10^(u)+(1-U)*10^(l)  # if instead want to sample space uniformly.

    return nothing
    
end

"""
    RPointSphereThetaPhi!()

Assigns the second (cos(theta)) and third (phi) elements of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta) and phi (phi normalised by pi). 
"""
function RPointSphereCosThetaPhi!(a::Vector{Float64}) 
    # Inputs a 4 element vector [p, cos(theta), phi/pi,theta/pi] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places)
    # phi points are normalised by pi

    u::Float64 = rand(Float64)
    v::Float64 = rand(Float64)

    a[2] = 2*v-1     # cos(theta) bound by [1,-1]
    a[3] = 2*u       # phi bound by [0,2) 
    a[4] = acos(a[2])/pi # theta bound by [0,1]

    return nothing
    
end

"""
    RPointSphereTheta!()

Assigns the second (cos(theta)) element of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta). 
"""
function RPointSphereCosTheta!(a::Vector{Float64}) 
    # Inputs a 3 element vector [p, cos(theta)] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta)
    # phi points are normalised by pi

    v::Float64 = rand(Float64)

    a[2] = 2*v-1     # cos(theta) bound by [1,-1]

    return nothing
    
end

"""
    betaVec!(βv,p1v,p2v)

Mutates the components of the centre of momentum velocity vector `βv` with components `[β,u,phi,γ,γβ] in terms of the incident state vectors `p1v` and `p1v`
"""
function betaVec!(βv::Vector{Float64},p1v,p2v,m1,m2)

    p1 = p1v[1]
    p2 = p2v[1]

    ct1 = p1v[2]
    ct2 = p2v[2]
    st1 = sqrt(1-ct1^2)
    st2 = sqrt(1-ct2^2)

    (sh1,ch1) = sincospi(p1v[3])
    (sh2,ch2) = sincospi(p2v[3])
    ch1h2 = cospi(p1v[3]-p2v[3])

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    E1::Float64 = Es1 + m1
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    E2::Float64 = Es2 + m2

    βv[1] = (1/(E1+E2)) * sqrt(p1^2 + p2^2 + 2*p1*p2*(ct1*ct2 + ch1h2*st1*st2))

    βv[2] = (1/(E1+E2)) * (p1*ct1 + p2*ct2) / βv[1] # cos(theta) bounded by [-1,1]
    x = p1*st1*ch1 + p2*st2*ch2
    y = p1*st1*sh1 + p2*st2*sh2
    βv[3] = mod(atan(y,x)/pi,2) # atan(y,x), note 1/(E1+E2) is common factor in y,x so ignored. phi bounded by [0,2]

    #gamma
    βv[4] = (E1+E2)
    βv[4] /= sqrt((m1+m2)^2 + 2*E1s*E2s+2*E1s*m2+2*E2s*m1-2*p1*p2*(ct1*ct2+ch1h2*st1*st2))

    # beta gamma
    βv[5] = sqrt(p1^2+p2^2+2*p1*p2*(ct1*ct2+ch1h2*st1*st2))
    βv[5] /= sqrt((m1+m2)^2 + 2*E1s*E2s+2*E1s*m2+2*E2s*m1-2*p1*p2*(ct1*ct2+ch1h2*st1*st2))
    
end

"""
    RPointSphereBoost!(βv,pvCOM)

Takes the random points on the sphere generated in the centre of momentum frame `pCv=[ct_deboosted,ct,h]` and boosts them back to the lab frame using `βv` to modify the lab frame vector `pLv`.
"""
function RPointSphereBoost!(pLv::Vector{Float64},βv::Vector{Float64},pCv::Vector{Float64})

    ctC = pCv[2]
    stC = sqrt(1-ctC^2)
    (shC,chC) = sincospi(pCv[3])

    ctβ = βv[2]
    stβ = sqrt(1-ctβ^2)
    (shβ,chβ) = sincospi(βv[3])

    pLv[2] = ctC*ctβ - stC*stβ*chC
    x = -stC*shC*shβ + chβ*(stC*chC*ctβ + ctC*stβ)
    y = stC*shC*chβ + shβ*(stC*chC*ctβ + ctC*stβ)

    pLv[3] = mod(atan(y,x)/pi,2)

    # deboost ctC to lab frame for probability calculation
    γ = βv[4]
    βγ = βv[5]
    #pCv[1] = (ctC+βv[1]) / (1+βv[1]*ctC)
    pCv[1] = (γ*ctC+βγ) / (γ+βγ*ctC)

end

"""
    pdfBoost(βv,pCv)

returns the probability of sampling a point given by `pCv` dependent on the boost `βv`.
"""
function pdfBoost(βv::Vector{Float64},ctC::Float64)

    # assumes pv in lab frame is valid and un-mirrored from the boosted frame
    β=βv[1]
    γ=βv[4]
    βγ=βv[5]

    #prob = (1-β^2)
    #prob /= 2*(1-β*ctC)^2
    
    prob = (γ^2-βγ^2)
    prob /= 2*(γ-βγ*ctC)^2
    
    return prob

end


"""
    WeightedFactors!(p1v,p2v,m1,m2,scale)

Returns the weighting rapidity `w` and the direction an angles `t` and `h` for rotations on a sphere. Weighting rapidity can be scaled by `scale`. Weight is dependant on the ratio of the energies of the two particles `p1v` and `p2v`. With angles `t` and `h` being the angles of the particle with the higher energy. 
"""
function WeightedFactors(p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,sBig::Float64,sSmol::Float64,scale::Float64) 

    E1::Float64 = sqrt(p1v[1]^2+m1^2)
    E2::Float64 = sqrt(p2v[1]^2+m2^2)
    ratio::Float64 = E1/E2
    gammaCOM::Float64 = (E1+E2)/sqrt(sBig+sSmol)

    if ratio > 1e2 # p1 has higher energy
        w3 = scale*log(ratio)
        t = p1v[4]
        h = p1v[3]
        # set w4 using (sBig+sSmol)/sBig i.e. how relativistic is the interaction is
        w4 = scale*log(gammaCOM)*4 # sqrt
    elseif ratio < 1e-2 # p2 has higher energy
        w4 = scale*log(1/ratio)
        t = p2v[4]
        h = p2v[3]
        # set w3 using (sBig+sSmol)/sBig i.e. how relativistic is the interaction is
        w3 = scale*log(gammaCOM)*4
    else # no weighting as uniform does well
        w3 = 0e0
        w4 = 0e0
        t = 0e0
        h = 0e0
    end

    return (w3,w4,t,h)
    
end

function WeightedFactors2(p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64,sBig::Float64,sSmol::Float64,scale::Float64) 

    p1 = p1v[1]
    p2 = p2v[1]

    (st1,ct1) = sincospi(p1v[4])
    (st2,ct2) = sincospi(p2v[4])
    (sh1,ch1) = sincospi(p1v[3])
    (sh2,ch2) = sincospi(p2v[3])
    ch1h2 = cospi(p1v[3]-p2v[3])

    s::Float64 = sBig + sSmol
    E1::Float64 = sqrt(p1v[1]^2+m1^2)
    E2::Float64 = sqrt(p2v[1]^2+m2^2)
    z::Float64 = sqrt(p1^2+p2^2+2*p1*p2*(ct1*ct2+ch1h2*st1*st2))
    x::Float64 = p1*st1*ch1 + p2*st2*ch2
    y::Float64 = p1*st1*sh1 + p2*st2*sh2

    # COM frame vector and angles
    γC::Float64 = (E1+E2)/sqrt(s)
    wC::Float64 = acosh(γC)
    tC::Float64 = acos((p1*ct1+p2*ct2)/z)/pi
    hC::Float64 = mod(atan(y,x)/pi,2)

    # outgoing COM frame momentum
    pC = InvariantFluxSmall(sSmol,m3,m4)/sqrt(s)

    # Lab Frame Angle limits
    # Due to conservation laws there is a limit on the lab frame scattering angle with respect to the COM frame direction. This angle is then mapped to a rapidity that "boosts" the lab frame angle samples such that 50% of samples lie within this angle limit.
    if m3 > m4
        w4Limit = 0e0
        tmp =  pC/(m3*sinh(wC))
        if tmp < 1e0 
            w3Limit = min(atanh(sqrt(1-tmp^2)),18.7e0) # for tmp < 1e-8 sqrt=0 due to float precision, 18.7e0 is maximum value of w3 to this precision
        else
            w3Limit = 0e0
        end
    elseif m3 < m4
        w3Limit = 0e0
        tmp = pC/(m4*sinh(wC))
        if tmp < 1e0
            w4Limit = min(atanh(sqrt(1-tmp^2)),18.7e0)
        else
            w4Limit = 0e0
        end
    else
        w3Limit = 0e0
        w4Limit = 0e0
    end

    wScale::Float64 = asinh(sqrt(sSmol/sBig)) # rough measure

    #w3::Float64 = min(w3Limit+wC+scale*wScale,18e0)
    #w4::Float64 = min(w4Limit+wC+scale*wScale,18e0)
    w3::Float64 = w3Limit+scale*(wC+wScale)
    #if p1v[1] < 2.0e0 && p2v[1] < 2.0e-10
    #    println("w3 = $w3")
    #    println("E1/E2 = $(log(E1/E2))")
    #end
    w4::Float64 = w4Limit+scale*(wC*wScale) 

    return (w3,w4,tC,hC)
    
end


"""
    WeightedFactorsEmission!(p1v,m1,scale)

Returns the weighting rapidity `w` and the direction an angles `t` and `h` for rotations on a sphere. Weighting rapidity can be scaled by `scale`. Weight is dependant on the energy of the emitting particle `p1v`. With angles `t` and `h` being the angles of the emitting particle. 
"""
function WeightedFactorsEmission(p1v::Vector{Float64},m1::Float64,scale::Float64) 

    E1::Float64 = sqrt(p1v[1]^2+m1^2)
    gamma::Float64 = E1/m1

    w = scale*acosh(gamma)
    t = p1v[4]
    h = p1v[3]

    return (w,t,h)
    
end

"""
    RPointSphereWeighted!()

Assigns randomly sampled angles on the sphere (cos(theta) and phi) weighed by a doppler boosting by rapidity `w`, returning the probability `prob` for such a sample and mutating the vector `a` with elements `[p, cos(theta), phi, theta]` with angles normalised by pi. 
"""
function RPointSphereWeighted!(a::Vector{Float64},w::Float64) 
    # Inputs a 5 element vector [p, cos(theta), phi, theta, prob] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html using the doppler boost formula as a weighting with a rapidity `w`. 
    # phi points are normalised by pi
    # prob is then the probability of sampling the random points on the sphere given the doppler boosting.

    v::Float64 = rand(Float64)
    #tB::Float64 = acos(2*v-1)/pi # "boosted" theta

    #costBdiv2sqr::Float64 = v
    #sintBdiv2sqr::Float64 = 1-v
    tantBdiv2::Float64 = sqrt((1-v)/v)


    h::Float64 = 2*rand(Float64) # phi points are normalised by pi
    #t::Float64 = 2*atan(exp(-w)*tanpi(tB/2))/pi  # "unboosted" theta divided by pi
    x::Float64 = exp(-w)*tantBdiv2
    t::Float64 = 2*atan(exp(-w)*tantBdiv2)/pi  # "unboosted" theta divided by pi
    
    #sintdiv2sqr::Float64 = x^2/(1+x^2)
    #costdiv2sqr::Float64 = 1/(1+x^2)

    ew::Float64 = exp(w)
    #prob::Float64 = (ew*sinpi(t/2)^2+cospi(t/2)^2/ew)^-2
    #prob::Float64 = (ew*sintdiv2sqr+costdiv2sqr/ew)^-2
    #prob::Float64 = (ew*cospi(tB/2)^2+sinpi(tB/2)^2/ew)^2
    #prob::Float64 = (ew*costBdiv2sqr+sintBdiv2sqr/ew)^2
    prob::Float64 = (ew*v+(1-v)/ew)^2

    a[2] = cospi(t)
    a[3] = h
    a[4] = t

#=     if w >= 7e0
        println("")
        println("w = $w")
        println("t = $t")
        println("tB = $tB")
        println("h = $h")
        println("prob = $prob")
        println("")
    end  =#

    return prob #/ (2*pi)
    
end

"""
    WeightedFactorsSync!()

Returns the weighting rapidity `w` and the direction an angles `t` and `h` for rotations on a sphere. Weighting rapidity can be scaled by `scale`. 
"""
function WeightedFactorsSync(pv::Vector{Float64},m::Float64,scale::Float64) 

    E1::Float64 = sqrt(pv[1]^2+m^2) # non dimensional units energy is gamma * m

    w = scale*acosh(E1/m) #*rand(Float64)
    t = pv[4]
    h = pv[3]

    return (w,t,h)
    
end