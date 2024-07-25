#= Functions for the S and T integation functions =#

"""
    TValue(p1v,p2v,sigma,mu1,mu2)

returns `Tval` with its Tval from MC integration based on initial momentum states `p1v` and `p2v` and cross section `sigma` based on particle selection. If initial state fails `sCheck`, i.e. cannot generate a physical output state, Tval is set to 0f0. 
Assumes f(x,p,μ)=constant over bin
"""
function TValue(p1v::Vector{Float32},p2v::Vector{Float32},sigma::Function,mu1::Float32,mu2::Float32,mu3::Float32,mu4::Float32)

    # define normalised masses
    m1 = mu1
    m2 = mu2
    m3 = mu3
    m4 = mu4

    # pre-defining terms for efficiency 
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct1::Float32 = p1v[2] #cospi(p1v[2])
    ct2::Float32 = p2v[2] #cospi(p2v[2]) 

    st1::Float32 = sqrt(1f0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float32 = sqrt(1f0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float32 = cospi(p1v[3]-p2v[3])

    m12 = m1^2
    m22 = m2^2
    
    #= This method leads to errors when p/m < 10^-6 due to Floating point precision, causing s<(m+m)^2 and InvFlux to retun a complex number error 
    E1 = sqrt(p1^2+m12)
    E2 = sqrt(p2^2+m22) 
    s = m12+m22 + 2*E1*E2 - 2*p1*p2*(ct1*ct2-ch1h2*st1*st2)
    =#

    #= will attempt to use the algebraic relation sqrt(1+x)-1 = x/(sqrt(1+x)+1) to split E up into a large and small part i.e. sqrt(p^2+m^2) = m + (p^2)/(sqrt(m^2+p^2)+m) = E 
    then label Es = (p^2)/(sqrt(m^2+p^2)+m) and rewrite 
    s = (m3+m4)^2 + 2(m3*Es4+m4*Es3+Es3*Es4 - p3p4(...)) =#

    Es1::Float32 = (p1^2)/(sqrt(m12+p1^2)+m1)
    Es1s::Float32 = Es1/p1
    Es2::Float32 = (p2^2)/(sqrt(m22+p2^2)+m2)
    Es2s::Float32 = Es2/p2

    sBig::Float32 = (m1+m2)^2
    sSmol::Float32 = 2*p1*p2*(- ct1*ct2-ch1h2*st1*st2 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    if sCheck(sSmol,sBig,m3,m4) # check if s value is valid for interaction
        E1::Float32 = Es1 + m1
        E2::Float32 = Es2 + m2
        
        # TSspe anisotropic absorption spectrum (to be integrated over d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
        Tval = (1/E1)*(1/E2)*(InvarientFluxSmall(sSmol,m1,m2))*sigma(sSmol,sBig)
        if (Tval==Inf32)
            error("ST1 Inf#"*string(sSmol)*"#"*string(sBig))
        end
    else # if not valid set T value to zero
        Tval = 0f0
    end

    return Tval

end

"""
    SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)

Returns `Sval` from MC integration based on initial momentum states `p1v` and `p2v` and final state `p3v` and differential cross section `dsigmadt` based on particle selection 12->34.  
Assumes f(x,p,μ)=constant over bin
"""
function SValue3(p3v::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32},dsigmadt::Function,mu1::Float32,mu2::Float32,mu3::Float32,mu4::Float32)

    # define normalise masses
    m1 = mu1
    m2 = mu2 
    m3 = mu3
    m4 = mu4 

    # pre-defining terms for efficiency 
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct1::Float32 = p1v[2] #cospi(p1v[2])
    ct2::Float32 = p2v[2] #cospi(p2v[2]) 

    st1::Float32 = sqrt(1f0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float32 = sqrt(1f0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float32 = cospi(p1v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2
    
    Es1::Float32 = m1 != 0f0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float32 = Es1/p1
    E1::Float32 = Es1 + m1
 
    Es2::Float32 = m2 != 0f0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float32 = Es2/p2
    E2::Float32 = Es2 + m2

    sBig::Float32 = (m1+m2)^2
    #sSmol::Float32 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    sSmol::Float32 = 2*p1*p2*(- ct1*ct2-ch1h2*st1*st2 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float32 = (1/E1)*(1/E2)*(2*InvarientFlux2Small(sSmol,m1,m2))

    p3::Float32 = p3v[1]
    ct3::Float32 = p3v[2] # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    st3::Float32 = sqrt(1f0-p3v[2]^2)
    ch3h1::Float32 = cospi(p3v[3]-p1v[3])
    ch3h2::Float32 = cospi(p3v[3]-p2v[3])
    Es3::Float32 = m3 != 0f0 ? (p3^2)/(sqrt(m32+p3^2)+m3) : p3
    Es3s::Float32 = Es3/p3

    # t = tBig + tSmol
    tBig::Float32 = (m3-m1)^2
    #tSmol::Float32 = -2*(m1*Es3 + m3*Es1 + Es3*Es1 - p3*p1*(ct3*ct1+ch3h1*st3*st1))
    tSmol::Float32 = -2*p3*p1*(- ct3*ct1 - ch3h1*st3*st1 + Es3s*Es1s + m1*Es3s/p1 + m3*Es1s/p3)
    # u = uBig + uSmol
    uBig::Float32 = (m2-m3)^2
    #uSmol::Float32 = m12+m22+m32+m42 - sBig - tBig - uBig - sSmol - tSmol  # this leads to Float32 issues better to calculate directly
    #uSmol::Float32 = -2*(m3*Es2 + m2*Es3 + Es2*Es3 - p2*p3*(ct2*ct3+ch3h2*st2*st3))
    uSmol::Float32 = -2*p2*p3*(-ct2*ct3 -ch3h2*st2*st3 + Es2s*Es3s + m3*Es2s/p3 + m2*Es3s/p2)

    deltacorrect::Float32 = (Es1*p3 - Es3*p1*(ct3*ct1+ch3h1*st3*st1) + Es2*p3 - Es3*p2*(ct3*ct2+ch3h2*st3*st2)) + (m1*p3 - m3*p1*(ct3*ct1+ch3h1*st3*st1) + m2*p3 - m3*p2*(ct3*ct2+ch3h2*st3*st2))
    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    Sval = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p3^2/(deltacorrect*sign(deltacorrect)))

    if (Sval==Inf || Sval == -Inf || Sval < 0f0)
        error("ST1 Inf#$deltacorrect#$sSmol#$sBig#$tSmol#$tBig#$uSmol#$uBig#$p3v#$p1v#$p2v")  
    end

    return Sval

end

"""
    SValue4(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)

Returns `Sval` from MC integration based on initial momentum states `p1v` and `p2v` and final state `p4v` and differential cross section `dsigmadt` based on particle selection 12->34.  
Assumes f(x,p,μ)=constant over bin
"""
function SValue4(p4v::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32},dsigmadt::Function,mu1::Float32,mu2::Float32,mu3::Float32,mu4::Float32)

    # define normalise masses
    m1 = mu1
    m2 = mu2 
    m3 = mu3
    m4 = mu4 

    # pre-defining terms for efficiency 
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct1::Float32 = p1v[2] #cospi(p1v[2])
    ct2::Float32 = p2v[2] #cospi(p2v[2]) 

    st1::Float32 = sqrt(1f0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float32 = sqrt(1f0-p2v[2]^2) #sinpi(p2v[2])

    ch1h2::Float32 = cospi(p1v[3]-p2v[3])

    m32 = m3^2
    m42 = m4^2
    m12 = m1^2
    m22 = m2^2

    Es1::Float32 = m1 != 0f0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    Es1s::Float32 = Es1/p1
    E1::Float32 = Es1 + m1
 
    Es2::Float32 = m2 != 0f0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    Es2s::Float32 = Es2/p2
    E2::Float32 = Es2 + m2

    sBig::Float32 = (m1+m2)^2
    #sSmol::Float32 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    sSmol::Float32 = 2*p1*p2*(- ct1*ct2 - ch1h2*st1*st2 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float32 = (1/E1)*(1/E2)*(2*InvarientFlux2Small(sSmol,m1,m2))

    p4::Float32 = p4v[1]
    ct4::Float32 = p4v[2] # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    st4::Float32 = sqrt(1f0-p4v[2]^2)
    ch4h1::Float32 = cospi(p4v[3]-p1v[3])
    ch4h2::Float32 = cospi(p4v[3]-p2v[3])
    Es4::Float32 = m4 != 0f0 ? (p4^2)/(sqrt(m42+p4^2)+m4) : p4
    Es4s::Float32 = Es4/p4
    E4::Float32 = Es4 + m4

    # u = uBig + uSmol
    uBig::Float32 = (m4-m1)^4
    #uSmol::Float32 = -2*(m1*Es4 + m4*Es1 + Es4*Es1 - p4*p1*(ct4*ct1+ch4h1*st4*st1))
    uSmol::Float32 = -2*p4*p1*(-ct4*ct1 - ch4h1*st4*st1 + Es4s*Es1s + m1*Es4s/p1 + m4*Es1s/p4)
    # t = tBig + tSmol
    tBig::Float32 = (m2-m4)^2
    #tSmol::Float32 = m12+m22+m32+m42 - sBig - uBig - tBig - sSmol - uSmol # Leads to Float32 issues, better to calculate directly 
    #tSmol::Float32 = -2*(m4*Es2 + m2*Es4 + Es2*Es4 - p2*p4*(ct2*ct4+ch4h2*st2*st4))
    tSmol::Float32 = -2*p2*p4*(-ct2*ct4 - ch4h2*st2*st4 + Es2s*Es4s + m4*Es2s/p4 + m2*Es4s/p2)

    deltacorrect::Float32 = (Es1*p4 - Es4*p1*(ct4*ct1+ch4h1*st4*st1) + Es2*p4 - Es4*p2*(ct4*ct2+ch4h2*st4*st2)) + (m1*p4 - m4*p1*(ct4*ct1+ch4h1*st4*st1) + m2*p4 - m4*p2*(ct4*ct2+ch4h2*st4*st2))
    # more Float accurate for when p1 and p2 have large order of magnitude difference as sum uses pairwise summation to reduce round of errors
    #sumTerms .= (Es1, -Es3prime*p1*(ct3*ct1+ch3h1*st3*st1), Es2, -Es3prime*p2*(ct3*ct2+ch3h2*st3*st2), m1, -(m3/p3)*p1*(ct3*ct1+ch3h1*st3*st1), m2, -(m3/p3)*p2*(ct3*ct2+ch3h2*st3*st2))
    #deltacorrect = p3*sum_oro(sumTerms)

    Sval = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p4^2/(deltacorrect*sign(deltacorrect)))

    if (Sval==Inf || Sval == -Inf || Sval < 0f0)
        error("ST1 Inf#$deltacorrect#$sSmol#$sBig#$tSmol#$tBig#$uSmol#$uBig#$p4v#$p1v#$p2v")  
    end

    return Sval

end

"""
    InvarientFlux(s,mu12,mu22)

returns the value of the invarient flux with 's' mandelstram variable and masses 'mass1' and 'mass2'
"""
function InvarientFlux(s::Float32,mu12::Float32,mu22::Float32)

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(s^2-2*s*(mu12+mu22)+(mu12-mu22)^2)/2

end

"""
    InvarientFluxSmall(sSmol,mu12,mu22)

returns the value of the invarient flux with smalled 's' mandelstram variable (sSmol = s - (m1+m2)^2)
"""
function InvarientFluxSmall(sSmol::Float32,mu1::Float32,mu2::Float32)
    # Better accuracy for small s

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(sSmol*(sSmol+4*mu1*mu2))/2

end

"""
    InvarientFlux2(s,mass12,mass22)

returns the value of the squared invarient flux with 's' mandelstram variable and masses 'mass1' and 'mass2'
"""
function InvarientFlux2(s::Float32,mu12::Float32,mu22::Float32)

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (s^2-2*s*(mu12+mu22)+(mu12-mu22)^2)/4

end

"""
    InvarientFluxSmall(sSmol,mass12,mass22)

returns the value of the squared invarient flux with smalled 's' mandelstram variable (sSmol = s - (m1+m2)^2)
"""
function InvarientFlux2Small(sSmol::Float32,mu1::Float32,mu2::Float32)
    # Better accuracy for small s

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (sSmol*(sSmol+4*mu1*mu2))/4

end




