    
"""
    Momentum3Value!(p3v,p3pv,p1v,p2v)

Takes set of random initial particle states 'p1v' and 'p2v' and random output states angles 'p3v[2:3]' and modifies outputs 'p3v' and 'p3pv' values with calculated output momentum and corrects angles if momentum is negative.
Function also returns a two bools 'p3_physical' and 'p3p_physical' indicating if p3 and p3p are physical momentum states given the inputs. 
Function also returns a Int 'NumStates' indicating the number of valid output states found.

# Examples
```julia-repl
julia> m1 = 1836.1528e0
julia> m2 = 1836.1528e0
julia> m3 = 1836.1528e0
julia> m4 = 1836.1528e0
julia> p1v = [1e0, 0.5e0, 1.8e0]
julia> p2v = [2e0, 0.2e0, 0.7e0]
julia> p3v = [0e0, 0.3e0, 0.7e0]
julia> p3pv = zeros(Float64,3)
julia> p3pv .= p3v
julia> Momentum3Value!(p3v,p3pv,p1v,p2v,m1,m2,m3,m4)
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
function Momentum3Value!(p3v::Vector{Float64},p3pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])
    st3::Float64,ct3::Float64 = sincospi(p3v[4])

    #ct3::Float64 = p3v[2] 
    #ct1::Float64 = p1v[2]
    #ct2::Float64 = p2v[2] 

    #st3::Float64 = sqrt(1e0-ct3^2)
    #st1::Float64 = sqrt(1e0-ct1^2)
    #st2::Float64 = sqrt(1e0-ct2^2)

    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    m32::Float64 = m3^2
    m42::Float64 = m4^2
    m12::Float64 = m1^2
    m22::Float64 = m2^2

    p12::Float64 = p1^2
    p22::Float64 = p2^2

    E1::Float64 = sqrt(m12+p12)
    E2::Float64 = sqrt(m22+p22)

    ct3ct1::Float64 = ct3*ct1 
    ct3ct2::Float64 = ct3*ct2
    ct1ct2::Float64 = ct1*ct2
    st3st1::Float64 = st3*st1
    st3st2::Float64 = st3*st2
    st1st2::Float64 = st1*st2

    # E1 = A1 + m1
    # E2 = A2 + m2
    A1::Float64 = p12/(E1+m1)
    A2::Float64 = p22/(E2+m2)

    #E1E2 = E1*E2

    #p1p2 = p1*p2

    # reset p3v values
    p3v[1] = 0e0 
    p3pv[1] = 0e0

    C3sqr::Float64 = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch3h1*st3st1+p2*ch3h2*st3st2))^2*(m32-m42+2*A2*m1+2*A1*A2+2*A1*m2+(m1+m2)^2-2*p1*p2*(ct1ct2+ch1h2*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch3h1*st3st1+p2*ch3h2*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch3h1*st3st1+p2*ch3h2*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*A2+2*A1*(-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch1h2*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch1h2*st1st2)) 

    C2::Float64 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch3h1*st3st1+p2*ch3h2*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch1h2*st1st2))

    C4::Float64 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch3h1*st3st1+p2*ch3h2*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch3h1*st3st1+p2*ch3h2*st3st2))

    # C3sqr == 0 was causing issues with SValue calculation often leading to deltacorrect = 0 so we are going to ignore this point and tread it as if p3 were complex.
    #=if C3sqr == 0 # only one state and p3 cannot equal zero

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
            if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
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

            if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2+m4
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

                if (p12/(E1+m1)+p22/(E2+m2)-p3p^2/(sqrt(m32+p3p^2)+m3)) > m3-m1-m2+m4
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
                if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2+m4
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

function Momentum3Value2!(p3v::Vector{Float64},p3pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # define identical states
    #identicalStates::Bool = false

    # pv should be [p,t,h]
    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]

    st1::Float64,ct1::Float64 = sincospi(p1v[4])
    st2::Float64,ct2::Float64 = sincospi(p2v[4])
    st3::Float64,ct3::Float64 = sincospi(p3v[4])

    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    m32::Float64 = m3^2
    m42::Float64 = m4^2
    m12::Float64 = m1^2
    m22::Float64 = m2^2

    p12::Float64 = p1^2
    p22::Float64 = p2^2

    E1::Float64 = sqrt(m12+p12)
    E2::Float64 = sqrt(m22+p22)

    ctheta12::Float64 = ct1*ct2 + ch1h2*st1*st2
    ctheta13::Float64 = ct3*ct1 + ch3h1*st3*st1
    ctheta23::Float64 = ct3*ct2 + ch3h2*st3*st2

    # E1 = A1 + m1
    # E2 = A2 + m2
    A1::Float64 = p12/(E1+m1)
    A2::Float64 = p22/(E2+m2)

    #E1E2 = E1*E2

    #p1p2 = p1*p2

    # reset p3v values
    p3v[1] = 0e0 
    p3pv[1] = 0e0

    C3sqr::Float64 = (E1+E2)^2*((m32-m42+m12+m22+2*E1*E2-2*p1*p2*ctheta12)^2-4*m32*((E1+E2)^2-(p1*ctheta13+p2*ctheta23)^2))

    C2::Float64 = -4*(p1*ctheta13+p2*ctheta23)*(m32-m42+m12+m22+2*E1*E2-2*p1*p2*ctheta12)

    C4::Float64 = -8*((E1+E2)^2-(p1*ctheta13+p2*ctheta23)^2)

    # C3sqr == 0 was causing issues with SValue calculation often leading to deltacorrect = 0 so we are going to ignore this point and tread it as if p3 were complex.
    #=if C3sqr == 0 # only one state and p3 cannot equal zero

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
            if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
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
                p3v[4] = 1-p3v[4]
            end

            if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2+m4
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
                    p3pv[4] = 1-p3pv[4]
                end

                if (p12/(E1+m1)+p22/(E2+m2)-p3p^2/(sqrt(m32+p3p^2)+m3)) > m3-m1-m2+m4
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
                    p3v[4] = 1-p3v[4]
                end
                if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2+m4
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
            p3v[4] = 1-p3v[4]
            #p3pv[2] *= -1
            #p3pv[3] = mod(p3pv[3]+1e0,2e0)
            end
        end

    end

    return p3_physical, p3p_physical, NumStates

end

function Momentum3Value_NoSignChange!(p3v::Vector{Float64},p3pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64)

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

    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])

    m32::Float64 = m3^2
    m42::Float64 = m4^2
    m12::Float64 = m1^2
    m22::Float64 = m2^2

    p12::Float64 = p1^2
    p22::Float64 = p2^2

    E1::Float64 = sqrt(m12+p12)
    E2::Float64 = sqrt(m22+p22)

    ct3ct1::Float64 = ct3*ct1 
    ct3ct2::Float64 = ct3*ct2
    ct1ct2::Float64 = ct1*ct2
    st3st1::Float64 = st3*st1
    st3st2::Float64 = st3*st2
    st1st2::Float64 = st1*st2

    A1::Float64 = p12/(E1+m1)
    A2::Float64 = p22/(E2+m2)

    #E1E2 = E1*E2

    #p1p2 = p1*p2

    # reset p3v values
    p3v[1] = 0e0 
    p3pv[1] = 0e0

    C3sqr::Float64 = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch3h1*st3st1+p2*ch3h2*st3st2))^2*(m32-m42+2*A2*m1+2*A1*A2+2*A1*m2+(m1+m2)^2-2*p1*p2*(ct1ct2+ch1h2*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch3h1*st3st1+p2*ch3h2*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch3h1*st3st1+p2*ch3h2*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*A2+2*A1*(-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch1h2*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch1h2*st1st2)) 

    C2::Float64 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch3h1*st3st1+p2*ch3h2*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch1h2*st1st2))

    C4::Float64 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch3h1*st3st1+p2*ch3h2*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch3h1*st3st1+p2*ch3h2*st3st2))

    # C3sqr == 0 was causing issues with SValue calculation often leading to deltacorrect = 0 so we are going to ignore this point and tread it as if p3 were complex.
    #=if C3sqr == 0 # only one state and p3 cannot equal zero

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
            if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2
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
                if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2+m4
                    p3_physical = true
                end
            else
                #p3v[1] = -p3
                #p3v[2] *= -1
                #p3v[3] = mod(p3v[3]+1e0,2e0)
                p3_physical = false
            end
        end

        if NumStates == 2
            p3p = (C2+C3)/C4
            if p3p == 0e0
                NumStates = 1
            else
                if p3p > 0e0
                    p3pv[1] = p3p
                    if (p12/(E1+m1)+p22/(E2+m2)-p3p^2/(sqrt(m32+p3p^2)+m3)) > m3-m1-m2+m4
                        p3p_physical = true
                    end
                else
                    #p3pv[1] = -p3p
                    #p3pv[2] *= -1
                    #p3pv[3] = mod(p3pv[3]+1e0,2e0)
                    p3_physical = false
                end
            end
        else # NumStates == 1
            p3 = (C2+C3)/C4
            if p3 == 0
                NumStates = 0
            else
                if p3 > 0e0
                    p3v[1] = p3
                    if (p12/(E1+m1)+p22/(E2+m2)-p3^2/(sqrt(m32+p3^2)+m3)) > m3-m1-m2+m4
                        p3_physical = true
                    end
                else
                    #p3v[1] = -p3
                    #p3v[2] *= -1
                    #p3v[3] = mod(p3v[3]+1e0,2e0)
                    p3_physical = false
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
        elseif p3Real < 0e0
            NumStates = 0
            #p3v[2] *= -1
            #p3v[3] = mod(p3v[3]+1e0,2e0)
        end

    end

    return p3_physical, p3p_physical, NumStates

end




"""
    pVector!(p4v,p3v,p1v,p2v)

Returns the p vector (in standard form [p,cos(theta),phi/pi]) of the non-sampled outgoing state, given the p1, p2 and p3 or p4 vectors using conservation of momentum.
"""
function pVector!(p4v::Vector{Float64},p3v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64})

    p1::Float64 = p1v[1]
    p2::Float64 = p2v[1]
    p3::Float64 = p3v[1]

    ct1::Float64 = p1v[2] 
    ct2::Float64 = p2v[2]  
    ct3::Float64 = p3v[2] 

    st1::Float64 = sqrt(1e0-p1v[2]^2)
    st2::Float64 = sqrt(1e0-p2v[2]^2) 
    st3::Float64 = sqrt(1e0-p3v[2]^2) 

    (sh1::Float64, ch1::Float64) = sincospi(p1v[3])
    (sh2::Float64, ch2::Float64) = sincospi(p2v[3])
    (sh3::Float64, ch3::Float64) = sincospi(p3v[3])
    ch1h2::Float64 = cospi(p1v[3]-p2v[3])
    ch3h1::Float64 = cospi(p3v[3]-p1v[3])
    ch3h2::Float64 = cospi(p3v[3]-p2v[3])

    a::Float64 = p1^2+p2^2+p3^2
    b::Float64 = 2*p1*p2*(ct1*ct2+ch1h2*st1*st2)
    c::Float64 = -2*p1*p3*(ct1*ct3+ch3h1*st1*st3)
    d::Float64 = -2*p2*p3*(ct2*ct3+ch3h2*st2*st3)
    p42::Float64 = a+b+c+d # making this one line seems to cause -ve values??
    #p42 = p1^1+p2^2+p3^2+2*p1*p2*(ct1*ct2+ch1h2*st1*st2)-2*p1*p3*(ct1*ct3+ch3h1*st1*st3)-2*p2*p3*(ct2*ct3+ch3h2*st2*st3)

    if p42 >= 0e0
        p4v[1] = sqrt(p42)
    else
        a = p1^2+p2^2+p3^2
        b = 2*p1*p2*(ct1*ct2+ch1h2*st1*st2)
        c = -2*p1*p3*(ct1*ct3+ch3h1*st1*st3)
        d = -2*p2*p3*(ct2*ct3+ch3h2*st2*st3)
        println("$a")
        println("$b")
        println("$c")
        println("$d")
        error("$p42")
    end

    p4v[2] = (p1*ct1+p2*ct2-p3*ct3)/p4v[1]

    if p4v[2] > 1e0
        ct4 = p4v[2]
        println("p4v[2] $ct4")
    end

    x = p1*st1*ch1+p2*st2*ch2-p3*st3*ch3 
    y = p1*st1*sh1+p2*st2*sh2-p3*st3*sh3

    p4v[3] = mod(atan(y,x)/pi,2)

    return nothing

end