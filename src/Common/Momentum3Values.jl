    
"""
    Momentum3Value!(p3v,p1v,p2v)

Takes set of random initial particle states 'p1v' and 'p2v' and random output states angles 'p3v[2:3,:]' and modifies output 'p3v[:,1]' values with calculated output momentum. Function also returns a bool 'identicalStates' as to whether the two output states are identical and bool 'testp3' 'testp3p' indicating if p3 and p3p are physical.

Requrires normalised masses (mu1,mu2,mu3,mu4) to be defined in advance in Init.jl.

# Examples
```julia-repl
julia> p1v = [1f0, 0.5f0, 1.8f0,]
julia> p2v = [2f0, 0.2f0, 0.7f0]
julia> p3v = [0f0 0f0; 0.3f0 0.3f0; 0.7f0 0.7f0]
julia> Momentum3Value!(p3v,p1v,p2v)
(false,true,true)
julia> p3v
 3-element Vector{Float32}:
 2.04505
 0.3
 0.7
julia> p3vp
 3-element Vector{Float32}
 0.691423
 -0.3
 1.7
```
"""
function Momentum3Value!(p3v::Array{Float32,2},p1v::Vector{Float32},p2v::Vector{Float32})

    # set normalised masses (defined in Init.jl)
    m1 = mu1
    m2 = mu2
    m3 = mu3
    m4 = mu4 

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct3::Float32 = p3v[2,1] #cospi(p3v[2,1]) # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    ct1::Float32 = p1v[2] #cospi(p1v[2])
    ct2::Float32 = p2v[2] #cospi(p2v[2]) 

    st3::Float32 = sqrt(1f0-ct3^2) #sinpi(p3v[2,1])
    st1::Float32 = sqrt(1f0-ct1^2) #sinpi(p1v[2])
    st2::Float32 = sqrt(1f0-ct2^2) #sinpi(p2v[2])

    ch1h3::Float32 = cospi(p3v[3,1]-p1v[3])
    ch1h4::Float32 = cospi(p3v[3,1]-p2v[3])
    ch3h4::Float32 = cospi(p1v[3]-p2v[3])

    m32::Float32 = m3^2
    m42::Float32 = m4^2
    m12::Float32 = m1^2
    m22::Float32 = m2^2

    p12::Float32 = p1^2
    p22::Float32 = p2^2

    sqm1p1::Float32 = sqrt(m12+p12)
    sqm2p2::Float32 = sqrt(m22+p22)

    ct3ct1::Float32 = ct3*ct1 
    ct3ct2::Float32 = ct3*ct2
    ct1ct2::Float32 = ct1*ct2
    st3st1::Float32 = st3*st1
    st3st2::Float32 = st3*st2
    st1st2::Float32 = st1*st2

    A1::Float32 = p12/(sqm1p1+m1)
    A2::Float32 = p22/(sqm2p2+m2)

    #sqm1p1sqm2p2 = sqm1p1*sqm2p2

    #p1p2 = p1*p2

    # reset p3v values
    p3v[1,1] = 0f0 
    p3v[1,2] = 0f0

    if ((C3sqr = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))^2*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*(A2-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))) > 0f0) # check for imaginary p3 state

        C2::Float32 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch3h4*st1st2))

        C3::Float32 = 4*sqrt(C3sqr)

        C4::Float32 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))

        val::Float32 = (C2-C3)/C4
        valp::Float32 = (C2+C3)/C4

        # either val or valp could be -ve or +ve and still allowed states, the statements below check each and adjust the angles appropriatly

        if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-val^2/(sqrt(m32+val^2)+m3)) > m3-m2-m1
            #assign for aligned case as valp +ve
            if val == 0f0
                testp3 = false
            elseif val > 0f0 # parallel
                testp3 = true
                p3v[1,1] = val 
            elseif valp < 0f0  # anti-parallel 
                testp3 = true
                p3v[1,1] = -val 
                p3v[2,1] *= -1 # mirrored in cos(theta) space is *-1. mod(1f0-p3v[2,1],1f0)     # theta bound by [0,1]
                p3v[3,1] = mod(p3v[3,1]+1f0,2f0)     # phi bound by [0,2) 
            else
                error("p3 state not accounted for"*string(val))
            end
        else
            testp3 = false
            #p3v[1,1] = 0f0
        end

        if (valp == val) # identical states C3 = 0 i.e. not to be counted
            NotIdenticalStates = false
        elseif (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-valp^2/(sqrt(m32+valp^2)+m3)) > m3-m2-m1 # non-identical but physical i.e. to be counted
            NotIdenticalStates = true
            #assign primed for aligned case as valp +ve
            if valp == 0f0
                testp3p = false
            elseif valp > 0f0 # parallel
                testp3p = true
                p3v[1,2] = valp 
            elseif valp < 0f0  # anti-parallel 
                testp3p = true
                p3v[1,2] = -valp 
                p3v[2,2] *= -1 # mod(1f0-p3v[2,2],1f0)     # theta bound by [0,1]
                p3v[3,2] = mod(p3v[3,2]+1f0,2f0)     # phi bound by [0,2) 
            else
                error("p3p state not accounted for:"*string(valp))
            end
        else # non-identical state but unphysical i.e. to be counted
            NotIdenticalStates = true
            testp3p = false
            #p3v[1,2] = 0f0
        end

    else # p3 is imagniary so unphysical but C3 != 0 so there are still two unphysical states
        testp3 = false
        testp3p = false
        NotIdenticalStates = true
        #p3v[1,1] = 0f0
        #p3v[1,2] = 0f0
    end

    return NotIdenticalStates, testp3, testp3p

end

function Momentum3Value2!(p3v::Vector{Float32},p3vp::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32})

    # set normalised masses (defined in Init.jl)
    m1 = mu1
    m2 = mu2
    m3 = mu3
    m4 = mu4 

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct3::Float32 = p3v[2] #cospi(p3v[2,1]) # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    ct1::Float32 = p1v[2] #cospi(p1v[2])
    ct2::Float32 = p2v[2] #cospi(p2v[2]) 

    st3::Float32 = sqrt(1f0-ct3^2) #sinpi(p3v[2,1])
    st1::Float32 = sqrt(1f0-ct1^2) #sinpi(p1v[2])
    st2::Float32 = sqrt(1f0-ct2^2) #sinpi(p2v[2])

    ch1h3::Float32 = cospi(p3v[3]-p1v[3])
    ch1h4::Float32 = cospi(p3v[3]-p2v[3])
    ch3h4::Float32 = cospi(p1v[3]-p2v[3])

    m32::Float32 = m3^2
    m42::Float32 = m4^2
    m12::Float32 = m1^2
    m22::Float32 = m2^2

    p12::Float32 = p1^2
    p22::Float32 = p2^2

    sqm1p1::Float32 = sqrt(m12+p12)
    sqm2p2::Float32 = sqrt(m22+p22)

    ct3ct1::Float32 = ct3*ct1 
    ct3ct2::Float32 = ct3*ct2
    ct1ct2::Float32 = ct1*ct2
    st3st1::Float32 = st3*st1
    st3st2::Float32 = st3*st2
    st1st2::Float32 = st1*st2

    A1::Float32 = p12/(sqm1p1+m1)
    A2::Float32 = p22/(sqm2p2+m2)

    #sqm1p1sqm2p2 = sqm1p1*sqm2p2

    #p1p2 = p1*p2

    # reset p3v values
    #p3v[1,1] = 0f0 
    #p3v[1,2] = 0f0

    C3sqr = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))^2*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*(A2-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2)) 

    C2::Float32 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch3h4*st1st2))

    C4::Float32 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))

    if C3sqr > 0f0 # p states are real and there are two of them

        NotIdenticalStates = true

        C3::Float32 = 4*sqrt(C3sqr)
        val::Float32 = (C2-C3)/C4
        valp::Float32 = (C2+C3)/C4

        # either val or valp could be -ve or +ve and still allowed states, the statements below check each and adjust the angles appropriatly

        # correct angles if val +/- and assign to p3v
        if val == 0f0      # zero momentum states are neglected
            zerop3 = true
        elseif val > 0f0 # parallel
            zerop3 = false
            p3v[1] = val 
        elseif val < 0f0  # anti-parallel 
            zerop3 = false
            p3v[1] = -val 
            p3v[2] *= -1 # mirrored in cos(theta) space is *-1. mod(1f0-p3v[2,1],1f0)     # theta bound by [0,1]
            p3v[3] = mod(p3v[3]+1f0,2f0)     # phi bound by [0,2) 
        else
            error("p3 state not accounted for"*string(val))
        end

        # check if val is physical
        if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-val^2/(sqrt(m32+val^2)+m3)) > m3-m2-m1 
            testp3 = true
        else
            testp3 = false
        end
        
        # correct angles if valp +/- and assign to p3vp
        if valp == 0f0              # zero momentum states are neglected
            zerop3p = true
        elseif valp > 0f0 # parallel
            zerop3p = false
            p3vp[1] = valp 
        elseif valp < 0f0  # anti-parallel 
            zerop3p = false
            p3vp[1] = -valp 
            p3vp[2] *= -1 # mod(1f0-p3v[2,2],1f0)     # theta bound by [0,1]
            p3vp[3] = mod(p3vp[3]+1f0,2f0)     # phi bound by [0,2) 
        else
            error("p3p state not accounted for:"*string(valp))
        end
        
        # check if valp physical
        if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-valp^2/(sqrt(m32+valp^2)+m3)) > m3-m2-m1            
            testp3p = true
        else
            testp3p = false
        end

    elseif C3sqr == 0f0 # only one state

        NotIdenticalStates = false
        testp3p = false
        zerop3p = false

        valIdentical::Float32 = C2/C4

        # correct angles if val +/- and assign to p3v
        if valIdentical == 0f0          # zero momentum states are neglected
            zerop3 = true
        elseif valIdentical > 0f0 # parallel
            zerop3 = false
            p3v[1] = valIdentical 
        elseif valIdentical < 0f0  # anti-parallel 
            zerop3 = false
            p3v[1] = -valIdentical 
            p3v[2] *= -1 # mirrored in cos(theta) space is *-1. mod(1f0-p3v[2,1],1f0)     # theta bound by [0,1]
            p3v[3] = mod(p3v[3]+1f0,2f0)     # phi bound by [0,2) 
        else
            error("p3 state not accounted for"*string(val))
        end

        # check if valp is physical
        if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-valIdentical^2/(sqrt(m32+valIdentical^2)+m3)) > m3-m2-m1 
            testp3 = true
        else
            testp3p = false
        end
        
    else # p3 is imagniary so one or two unphysical states

        NotIdenticalStates = true
        testp3 = false
        testp3p = false

        # need to correctly assign t3 values depending on if real part of p3 is +/-
        valReal = C2/C4
        # if valReal + then no change, if -ve then need to re-assign costheta and phi for both states
        if valReal == 0f0 # real part is zero but imaginary part is not so two states
            zerop3 = true
            zerop3p = true
        elseif valReal < 0f0
            zerop3 = false
            zerop3p = false
            p3v[2] *= -1
            p3vp[2] *= -1
            p3v[3] = mod(p3v[3]+1f0,2f0)
            p3vp[3] = mod(p3vp[3]+1f0,2f0)
        else
            zerop3 = false
            zerop3p = false
        end
        
    end

    return NotIdenticalStates, testp3, testp3p, zerop3, zerop3p

end

function Momentum3Value3!(p3v::Vector{Float32},p3pv::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32})

    # set normalised masses (defined in Init.jl)
    m1 = mu1
    m2 = mu2
    m3 = mu3
    m4 = mu4 

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct3::Float32 = p3v[2] 
    ct1::Float32 = p1v[2]
    ct2::Float32 = p2v[2] 

    st3::Float32 = sqrt(1f0-ct3^2)
    st1::Float32 = sqrt(1f0-ct1^2)
    st2::Float32 = sqrt(1f0-ct2^2)

    ch1h3::Float32 = cospi(p3v[3]-p1v[3])
    ch1h4::Float32 = cospi(p3v[3]-p2v[3])
    ch3h4::Float32 = cospi(p1v[3]-p2v[3])

    m32::Float32 = m3^2
    m42::Float32 = m4^2
    m12::Float32 = m1^2
    m22::Float32 = m2^2

    p12::Float32 = p1^2
    p22::Float32 = p2^2

    sqm1p1::Float32 = sqrt(m12+p12)
    sqm2p2::Float32 = sqrt(m22+p22)

    ct3ct1::Float32 = ct3*ct1 
    ct3ct2::Float32 = ct3*ct2
    ct1ct2::Float32 = ct1*ct2
    st3st1::Float32 = st3*st1
    st3st2::Float32 = st3*st2
    st1st2::Float32 = st1*st2

    A1::Float32 = p12/(sqm1p1+m1)
    A2::Float32 = p22/(sqm2p2+m2)

    #sqm1p1sqm2p2 = sqm1p1*sqm2p2

    #p1p2 = p1*p2

    # reset p3v values
    #p3v[1,1] = 0f0 
    #p3v[1,2] = 0f0

    C3sqr::Float32 = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))^2*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*(A2-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2)) 

    C2::Float32 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch3h4*st1st2))

    C4::Float32 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))

    if C3sqr == 0f0 # only one state
        NotIdentical = false
        p3p_physical = false

        p3 = C2/C4
        if p3 >= 0f0
            p3v[1] = p3
        elseif p3 < 0f0
            p3v[1] = -p3
            p3v[2] *= -1
            p3v[3] = mod(p3v[3]+1f0,2f0)
        end
        if p3 != 0f0 && (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m2-m1
            p3_physical = true
        else
            p3_physical = false
        end

    elseif C3sqr > 0f0
        NotIdentical = true
        C3 = sqrt(C3sqr)
        p3 = (C2-C3)/C4
        p3p = (C2+C3)/C4
        if p3 >= 0f0
            p3v[1] = p3
        else
            p3v[1] = -p3
            p3v[2] *= -1
            p3v[3] = mod(p3pv[3]+1f0,2f0)
        end
        if p3p >= 0f0
            p3pv[1] = p3p
        else
            p3pv[1] = -p3p
            p3pv[2] *= -1
            p3pv[3] = mod(p3pv[3]+1f0,2f0)
        end
        if p3 != 0f0 && (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m2-m1
            p3_physical = true
        else
            p3_physical = false
        end
        if p3p != 0f0 && (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3p^2/(sqrt(m32+p3p^2)+m3)) > m3-m2-m1
            p3p_physical = true
        else
            p3p_physical = false
        end

    else # imaginary

        NotIdentical = true
        p3p_physical = false
        p3_physical = false

        p3Real = C2/C4
        if p3Real < 0f0
            p3v[2] *= -1
            p3v[3] = mod(p3pv[3]+1f0,2f0)
        end

    end
    
    #=
    C3::ComplexF32 = 4*sqrt(complex(C3sqr))

    p3Complex::ComplexF32 = (C2-C3)/C4
    p3pComplex::ComplexF32 = (C2+C3)/C4

    p3::Float32 = real(p3Complex)
    p3sign::Bool = signbit(p3) 

    # direction adjustments
    if p3sign # true if -ve
        p3v[1] = -p3
        p3v[2] *= -1
        p3v[3] = mod(p3v[3]+1f0,2f0)
    else
        p3v[1] = p3
    end
    # physical checks
    if isreal(p3Complex) && (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-real(p3)^2/(sqrt(m32+real(p3)^2)+m3)) > m3-m2-m1 
        p3_physical = true
    else
        p3_physical = false
    end
    # are they identical
    if C3sqr != 0f0
        NotIdentical = true
        p3p::Float32 = real(p3pComplex)
        p3psign::Bool = signbit(p3p)
        if p3psign # true if -ve
            p3pv[1] = -p3p
            p3pv[2] *= -1
            p3pv[3] = mod(p3pv[3]+1f0,2f0)
        else
            p3pv[1] = p3p
        end
        # physical check
        if isreal(p3pComplex) && (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-real(p3p)^2/(sqrt(m32+real(p3p)^2)+m3)) > m3-m2-m1 
            p3p_physical = true
        else
            p3p_physical = false
        end
    else
        NotIdentical = false
        p3p_physical = false
    end
    =#

    return p3_physical, p3p_physical, NotIdentical

end


#= Testing
p1v = [1f0, 0.5f0, 1.8f0,]
p2v = [2f0, 0.2f0, 0.7f0]
p3v = [0f0 0f0; 0.3f0 0.3f0; 0.7f0 0.7f0]

Momentum3Value!(p3v,p1v,p2v)
p3v
=#