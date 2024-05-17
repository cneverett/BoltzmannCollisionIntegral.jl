#= Functions for the S and T integation functions =#

function TValue!(ST::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32},m1::Float32,m2::Float32)

    # returns one t value

    # pre-defining terms for efficiency 
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct1::Float32 = cospi(p1v[2])
    ct2::Float32 = cospi(p2v[2]) 

    st1::Float32 = sinpi(p1v[2])
    st2::Float32 = sinpi(p2v[2])

    ch1h2::Float32 = cospi(p1v[3]-p2v[3])

    m12 = m1^2
    m22 = m2^2
    
    #= This method leads to errors when p/m < 10^-6 due to floating point precision, causing s<(m+m)^2 and InvFlux to retun a complex number error 
    E1 = sqrt(p1^2+m12)
    E2 = sqrt(p2^2+m22) 
    s = m12+m22 + 2*E1*E2 - 2*p1*p2*(ct1*ct2-ch1h2*st1*st2)
    =#

    #= will attempt to use the algebraic relation sqrt(1+x)-1 = x/(sqrt(1+x)+1) to split E up into a large and small part i.e. sqrt(p^2+m^2) = m + (p^2)/(sqrt(m^2+p^2)+m) = E 
    then label Es = (p^2)/(sqrt(m^2+p^2)+m) and rewrite 
    s = (m3+m4)^2 + 2(m3*Es4+m4*Es3+Es3*Es4 - p3p4(...)) =#

    Es1::Float32 = (p1^2)/(sqrt(m12+p1^2)+m1)
    Es2::Float32 = (p2^2)/(sqrt(m22+p2^2)+m2)

    sBig::Float32 = (m1+m2)^2
    sSmol::Float32 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2-ch1h2*st1*st2))

    #s = sBig + sSmol
    
    E1::Float32 = Es1 + m1
    E2::Float32 = Es2 + m2
    
    #println(s)
    #println(t)

    #= if InvFlux2(s,m32,m42)==0f0
        println("s:"*string(s)*"ss:"*string(sSmol)*"#sb:"*string(sBig)*"#if2s:"*string(InvFlux2Small(sSmol,m3,m4)))
    end =#

    # TSspe anisotropic absorption spectrum (to be integrated over d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    ST[3] = (1/E1)*(1/E2)*(InvFluxSmall(sSmol,m1,m2))*sigma(sSmol,sBig)
    #println(sBig)
    #println(sSmol)
    #println(sigma(sSmol,sBig))
    if (ST[3]==Inf32)
        println(string)
        ST[3] = 0f0
    end

    #println(string(ST[1],ST[2]))

    return nothing

end

function TValuewithTest!(ST::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32},m1::Float32,m2::Float32)

    # returns one t value

    # pre-defining terms for efficiency 
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct1::Float32 = cospi(p1v[2])
    ct2::Float32 = cospi(p2v[2]) 

    st1::Float32 = sinpi(p1v[2])
    st2::Float32 = sinpi(p2v[2])

    ch1h2::Float32 = cospi(p1v[3]-p2v[3])

    m12 = m1^2
    m22 = m2^2
    
    #= This method leads to errors when p/m < 10^-6 due to floating point precision, causing s<(m+m)^2 and InvFlux to retun a complex number error 
    E1 = sqrt(p1^2+m12)
    E2 = sqrt(p2^2+m22) 
    s = m12+m22 + 2*E1*E2 - 2*p1*p2*(ct1*ct2-ch1h2*st1*st2)
    =#

    #= will attempt to use the algebraic relation sqrt(1+x)-1 = x/(sqrt(1+x)+1) to split E up into a large and small part i.e. sqrt(p^2+m^2) = m + (p^2)/(sqrt(m^2+p^2)+m) = E 
    then label Es = (p^2)/(sqrt(m^2+p^2)+m) and rewrite 
    s = (m3+m4)^2 + 2(m3*Es4+m4*Es3+Es3*Es4 - p3p4(...)) =#

    Es1::Float32 = (p1^2)/(sqrt(m12+p1^2)+m1)
    Es2::Float32 = (p2^2)/(sqrt(m22+p2^2)+m2)

    sBig::Float32 = (m1+m2)^2
    sSmol::Float32 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))

    if sCheck(sSmol,sBig) # check if s value is valid for interaction
        E1::Float32 = Es1 + m1
        E2::Float32 = Es2 + m2
        
        #println(s)
        #println(t)

        #= if InvFlux2(s,m32,m42)==0f0
            println("s:"*string(s)*"ss:"*string(sSmol)*"#sb:"*string(sBig)*"#if2s:"*string(InvFlux2Small(sSmol,m3,m4)))
        end =#

        # TSspe anisotropic absorption spectrum (to be integrated over d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
        ST[3] = (1/E1)*(1/E2)*(InvFluxSmall(sSmol,m1,m2))*sigma(sSmol,sBig)
        #println(sBig)
        #println(sSmol)
        #println(sigma(sSmol,sBig))
        if (ST[3]==Inf32)
            println(sSmol,sBig)
            ST[3] = 0f0
        end
    else # if not valid set T value to zero
        ST[3] = 0f0
    end

    return nothing

end

function SValueWithTests!(ST::Vector{Float32},p3v::Array{Float32},p1v::Vector{Float32},p2v::Vector{Float32},m1::Float32,m2::Float32,m3::Float32,testp3::Bool,testp3p::Bool)

    # returns two s values

    # pre-defining terms for efficiency 
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct1::Float32 = cospi(p1v[2])
    ct2::Float32 = cospi(p2v[2]) 

    st1::Float32 = sinpi(p1v[2])
    st2::Float32 = sinpi(p2v[2])

    ch1h2::Float32 = cospi(p1v[3]-p2v[3])

    m32 = m3^2
    m12 = m1^2
    m22 = m2^2
    
    #= This method leads to errors when p/m < 10^-6 due to floating point precision, causing s<(m+m)^2 and InvFlux to retun a complex number error 
    E1 = sqrt(p1^2+m12)
    E3 = sqrt(p3^2+m32)
    E4 = sqrt(p4^2+m42) 
    s = m32+m42 + 2*E3*E4 - 2*p3*p4*(ct3*ct4-ch3h4*st3*st4)
    t = m32+m12 - 2*E1*E3 + 2*p1*p3*(ct1*ct3+ch1h3*st1*st3)
    =#

    #= will attempt to use the algebraic relation sqrt(1+x)-1 = x/(sqrt(1+x)+1) to split E up into a large and small part i.e. sqrt(p^2+m^2) = m + (p^2)/(sqrt(m^2+p^2)+m) = E 
    then label Es = (p^2)/(sqrt(m^2+p^2)+m) and rewrite 
    s = (m3+m4)^2 + 2(m3*Es4+m4*Es3+Es3*Es4 - p3p4(...))
    t = (m1-m3)^2 - 2(m3*Es1+m1*Es3+Es1*Es3 - p1p3(...)) =#

    Es1::Float32 = (p1^2)/(sqrt(m12+p1^2)+m1)
    Es2::Float32 = (p2^2)/(sqrt(m22+p2^2)+m2)

    sBig::Float32 = (m1+m2)^2
    sSmol::Float32 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))

    #s = sBig + sSmol
    
    E1::Float32 = Es1 + m1
    E2::Float32 = Es2 + m2
    
    #println(s)
    #println(t)

    #= if InvFlux2(s,m32,m42)==0f0
        println("s:"*string(s)*"ss:"*string(sSmol)*"#sb:"*string(sBig)*"#if2s:"*string(InvFlux2Small(sSmol,m3,m4)))
    end =#

    # Sspe anisotropic emission spectrum (to be integrated over d^2p1d^3p3d^3p4). See obsidian note on discrete anisotropic kinetic equation
    val::Float32 = (1/E1)*(1/E2)*(2*InvFlux2Small(sSmol,m1,m2))
    #println(sBig)
    #println(sSmol)
    #println(sigma(sSmol,sBig))

    if testp3   # t value for p3
        p3::Float32 = p3v[1,1]
        ct3::Float32 = cospi(p3v[2,1]) # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
        st3::Float32 = sinpi(p3v[2,1])
        ch3h1::Float32 = cospi(p3v[3,1]-p1v[3])
        ch3h2::Float32 = cospi(p3v[3,1]-p2v[3])
        Es3::Float32 = (p3^2)/(sqrt(m32+p3^2)+m3)

        # t = tBig + tSmol
        tBig::Float32 = (m3-m1)^2
        tSmol::Float32 = -2*(m1*Es3 + m3*Es1 + Es3*Es1 - p3*p1*(ct3*ct1+ch3h1*st3*st1))
        # u = uBig + uSmol
        uBig::Float32 = (m2-m3)^2
        uSmol::Float32 = -2*(m3*Es2 + m2*Es3 + Es2*Es3 - p2*p3*(ct2*ct3+ch3h2*st2*st3))

        deltacorrect::Float32 = (Es1*p3 - Es3*p1*(ct3*ct1+ch3h1*st3*st1) + Es2*p3 - Es3*p2*(ct3*ct2+ch3h2*st3*st2)) + (m1*p3 - m3*p1*(ct3*ct1+ch3h1*st3*st1) + m2*p3 - m3*p2*(ct3*ct2+ch3h2*st3*st2))

        #if (uCheck(uSmol,uBig)&&tCheck(tSmol,tBig)#= &&stuCheck(sSmol,sBig,tSmol,tBig,uSmol,uBig) =#) == false
        #    error("stu#"*string(tSmol)*"#"*string(tBig)*"#"*string(uSmol)*"#"*string(uBig)*"#"*string(sSmol)*"#"*string(sSmol*(tSmol+tBig)*(uSmol+uBig)-sSmol)*"#"*string(sSmol+(tSmol+uSmol)+(tBig+uBig)))
        #end

        ST[1] = dsigmadt(sSmol,sBig,tSmol,tBig,uSmol,uBig)*val*(p3^2/(deltacorrect*sign(deltacorrect)))
        if (ST[1]==Inf || ST[1] == -Inf)
            error("ST1 Inf#"*string(deltacorrect)*"#"*string(tSmol)*"#"*string(tBig)*"#"*string(sSmol)*"#"*string(sBig))  
        end

        #if dsigmadt(s,t) == 0f0
        #    println("s,t:"*string(s)*","*string(t))
        #end

        #println(string(tSmol+tBig)*"#"*string(uSmol+uBig)*"#"*string(sSmol))
        #println(sSmol+tSmol+tBig+uSmol+uBig)

    end

    if testp3p # t value of p3p
        p3p::Float32 = p3v[1,2]
        ct3p::Float32 = cospi(p3v[2,2]) # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
        st3p::Float32 = sinpi(p3v[2,2])
        ch3h1p::Float32 = cospi(p3v[3,2]-p1v[3])
        ch3h2p::Float32 = cospi(p3v[3,2]-p2v[3])
        Es3p::Float32 = (p3p^2)/(sqrt(m32+p3p^2)+m3)

        tBigp::Float32 = (m3-m1)^2
        tSmolp::Float32 = -2*(m1*Es3p + m3*Es1 + Es3p*Es1 - p3p*p1*(ct3p*ct1+ch3h1p*st3p*st1))
        uBigp::Float32 = (m2-m3)^2
        uSmolp::Float32 = -2*(m3*Es2 + m2*Es3p + Es2*Es3p - p2*p3p*(ct2*ct3p+ch3h2p*st2*st3p))
        deltacorrectp::Float32 = (Es1*p3p - Es3p*p1*(ct3p*ct1+ch3h1p*st3p*st1) + Es2*p3p - Es3p*p2*(ct3p*ct2+ch3h2p*st3p*st2)) + (m1*p3p - m3*p1*(ct3p*ct1+ch3h1p*st3p*st1) + m2*p3p - m3*p2*(ct3p*ct2+ch3h2p*st3p*st2))

        #if (uCheck(uSmolp,uBigp)&&tCheck(tSmolp,tBigp)#= &&stuCheck(sSmol,sBig,tSmolp,tBigp,uSmolp,uBigp) =#) == false
        #    error("stup#"*string(tSmolp+tBigp)*"#"*string(uSmolp+uBigp)*"#"*string(sSmol)*"#"*string(sSmol*(tSmolp+tBigp)*(uSmolp+uBigp)-sSmol)*"#"*string(sSmol+tSmolp+uSmolp+uBigp+tBigp))
        #end

        ST[2] = dsigmadt(sSmol,sBig,tSmolp,tBigp,uSmolp,uBigp)*val*(p3p^2/(deltacorrectp*sign(deltacorrectp)))
        if (ST[2]==Inf || ST[2] == -Inf)
            error("ST2 Inf#"*string(deltacorrectp)*"#"*string(tSmolp)*"#"*string(tBigp))
        end

        #if dsigmadt(s,tp) == 0f0
        #    println("s,t:"*string(s)*","*string(tp))
        #end

        #println(sSmol+tSmolp+tBigp+uSmolp+uBigp)

    end

    #println(string(ST[1],ST[2]))

    return nothing

end

function InvFlux(s::Float32,mass12::Float32,mass22::Float32)

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(s^2-2*s*(mass12+mass22)+(mass12-mass22)^2)/2


end

function InvFluxSmall(sSmol::Float32,mass1::Float32,mass2::Float32)
    # Better accuracy for small s

    # sqrt(lambda(s,m1^2,m2^2))/2 = sqrt(s)|p*|
    return sqrt(sSmol*(sSmol+4*mass1*mass2))/2

end

function InvFlux2(s::Float32,mass12::Float32,mass22::Float32)

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (s^2-2*s*(mass12+mass22)+(mass12-mass22)^2)/4

end

function InvFlux2Small(sSmol::Float32,mass1::Float32,mass2::Float32)
    # Better accuracy for small s

    # lambda(s,m1^2,m2^2)/4 = s|p*|^2
    return (sSmol*(sSmol+4*mass1*mass2))/4

end


