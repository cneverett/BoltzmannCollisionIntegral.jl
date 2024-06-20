    
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
 julia p3vp
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

    st3::Float32 = sqrt(1f0-p3v[2,1]^2) #sinpi(p3v[2,1])
    st1::Float32 = sqrt(1f0-p1v[2]^2) #sinpi(p1v[2])
    st2::Float32 = sqrt(1f0-p2v[2]^2) #sinpi(p2v[2])

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

        if (m1+m2-m3+p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-val^2/(sqrt(m32+val^2)+m3))>0
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
        elseif ((m1+m2-m3+p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-valp^2/(sqrt(m32+valp^2)+m3))>0) # non-identical but physical i.e. to be counted
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


#= Testing
p1v = [1f0, 0.5f0, 1.8f0,]
p2v = [2f0, 0.2f0, 0.7f0]
p3v = [0f0 0f0; 0.3f0 0.3f0; 0.7f0 0.7f0]

Momentum3Value!(p3v,p1v,p2v)
p3v
=#