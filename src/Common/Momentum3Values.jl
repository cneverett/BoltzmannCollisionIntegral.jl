    
"""
    Momentum3Value!(p3v,p3pv,p1v,p2v)

Takes set of random initial particle states 'p1v' and 'p2v' and random output states angles 'p3v[2:3]' and modifies outputs 'p3v' and 'p3pv' values with calculated output momentum and corrects angles if momentum is negative.
Function also returns a two bools 'p3_physical' and 'p3p_physical' indicating if p3 and p3p are physical momentum states given the inputs. 
Function also returns a Int 'NumStates' indicating the number of valid output states found.

Requrires normalised masses (mu1,mu2,mu3,mu4) to be defined in advance in Init.jl as const.

# Examples
```julia-repl
julia> mu1 = 1836.1528e0
julia> mu2 = 1836.1528e0
julia> mu3 = 1836.1528e0
julia> mu4 = 1836.1528e0
julia> p1v = [1e0, 0.5e0, 1.8e0]
julia> p2v = [2e0, 0.2e0, 0.7e0]
julia> p3v = [0e0, 0.3e0, 0.7e0]
julia> p3pv = zeros(Float64,3)
julia> p3pv .= p3v
julia> Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)
(true,true,2)
julia> p3v
 3-element Vector{Float64}:
 2.04505
 0.3
 0.7
julia> p3pv
 3-element Vector{Float64}
 0.691423
 -0.3
 1.7
```
"""
function Momentum3Value!(p3v::Vector{Float64},p3pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)

    # set normalised masses (defined in Init.jl)
    m1 = mu1
    m2 = mu2
    m3 = mu3
    m4 = mu4 

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    ct3::Float64 = p3v[2] 
    ct1::Float64 = p1v[2]
    ct2::Float64 = p2v[2] 

    st3::Float64 = sqrt(1e0-ct3^2)
    st1::Float64 = sqrt(1e0-ct1^2)
    st2::Float64 = sqrt(1e0-ct2^2)

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

    A1::Float64 = p12/(sqm1p1+m1)
    A2::Float64 = p22/(sqm2p2+m2)

    #sqm1p1sqm2p2 = sqm1p1*sqm2p2

    #p1p2 = p1*p2

    # reset p3v values
    p3v[1] = 0e0 
    p3pv[1] = 0e0

    C3sqr::Float64 = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))^2*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*(A2-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2)) 

    C2::Float64 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch3h4*st1st2))

    C4::Float64 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))

    # C3sqr == 0 was causing issues with SValue calculation often leading to deltacorrect = 0 so we are going to ignore this point and tread it as if p3 were complex.
    #=if C3sqr == 0 # only one state and p3 cannont equal zero

        NumStates = 1
        p3_physical = false
        p3p_physical = false

        p3 = C2/C4
        if p3 == 0e0
            NumStates = 0
        else
            if p3 > 0e0
                p3v[1] = p3
            else
                p3v[1] = -p3
                p3v[2] *= -1
                p3v[3] = mod(p3v[3]+1e0,2e0)
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

        if p3 == 0e0
            NumStates = 1
        else
            if p3 > 0e0
                p3v[1] = p3
            else
                p3v[1] = -p3
                p3v[2] *= -1
                p3v[3] = mod(p3v[3]+1e0,2e0)
            end

            if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
                p3_physical = true
            end
        end

        if NumStates == 2
            p3p = (C2+C3)/C4
            if p3p == 0e0
                NumStates = 1
            else
                if p3p > 0e0
                    p3pv[1] = p3p
                    
                else
                    p3pv[1] = -p3p
                    p3pv[2] *= -1
                    p3pv[3] = mod(p3pv[3]+1e0,2e0)
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
                if p3 > 0e0
                    p3v[1] = p3
                else
                    p3v[1] = -p3
                    p3v[2] *= -1
                    p3v[3] = mod(p3v[3]+1e0,2e0)
                end
                if (p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
                    p3_physical = true
                end
            end    

        end

    else # imaginary C3sqr < 0e0

        NumStates = 1 # two states but both in same bin so same as one
        p3p_physical = false
        p3_physical = false

        p3Real = C2/C4
        if p3Real == 0e0
            NumStates = 0
        else
            if p3Real < 0e0
            p3v[2] *= -1
            p3v[3] = mod(p3v[3]+1e0,2e0)
            #p3pv[2] *= -1
            #p3pv[3] = mod(p3pv[3]+1e0,2e0)
            end
        end

    end

    return p3_physical, p3p_physical, NumStates

end


"""
    p4Vector!(p4v,p3v,p1v,p2v)

Returns the p4 vector (in standard form [p,cos(theta),phi/pi]) given the p1, p2 and p3 vectors using conservation of momentum.
"""
function p4Vector!(p4v::Vector{Float64},p3v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64})

    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    p3::Float64 = p3v[1]

    ct1::Float64 = p1v[2] 
    ct2::Float64 = p2v[2]  
    ct3::Float64 = p3v[2] 

    st1::Float64 = sqrt(1e0-p1v[2]^2)
    st2::Float64 = sqrt(1e0-p2v[2]^2) 
    st3::Float64 = sqrt(1e0-p3v[2]^2) 

    ch1::Float64 = cospi(p1v[3])
    ch2::Float64 = cospi(p2v[3])
    ch3::Float64 = cospi(p3v[3])

    sh1::Float64 = sqrt(1e0-ch1^2) 
    sh2::Float64 = sqrt(1e0-ch2^2) 
    sh3::Float64 = sqrt(1e0-ch3^2) 

    p3xyz = [p3*st3*ch3,p3*st3*sh3,p3*ct3]
    p1xyz = [p1*st1*ch1,p1*st1*sh1,p1*ct1]
    p2xyz = [p2*st2*ch2,p2*st2*sh2,p2*ct2]

    p4xyz = p1xyz + p2xyz - p3xyz

    p4v[1] = sqrt(p4xyz[1]^2+p4xyz[2]^2+p4xyz[3]^2)
    p4v[2] = p4xyz[3]/p4v[1]
    p4v[3] = atan(p4xyz[2],p4xyz[1])/pi

    return nothing

end