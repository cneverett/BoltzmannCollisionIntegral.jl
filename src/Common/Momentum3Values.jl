    
"""
    Momentum3Value!(p3v,p3pv,p1v,p2v)

Takes set of random initial particle states 'p1v' and 'p2v' and random output states angles 'p3v[2:3]' and modifies outputs 'p3v' and 'p3pv' values with calculated output momentum and corrects angles if momentum is negative.
Function also returns a two bools 'p3_physical' and 'p3p_physical' indicating if p3 and p3p are physical momentum states given the inputs. 
Function also returns a Int 'NumStates' indicating the number of valid output states found.

Requrires normalised masses (mu1,mu2,mu3,mu4) to be defined in advance in Init.jl as const.

# Examples
```julia-repl
julia> p1v = [1f0, 0.5f0, 1.8f0]
julia> p2v = [2f0, 0.2f0, 0.7f0]
julia> p3v = [0f0, 0.3f0, 0.7f0]
julia> p3pv = zeros(Float32,3)
julia> p3pv .= p3v
julia> Momentum3Value!(p3v,p3pv,p1v,p2v)
(true,true,2)
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
function Momentum3Value!(p3v::Vector{Float32},p3pv::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32})

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
    p3v[1] = 0f0 
    p3pv[1] = 0f0

    C3sqr::Float32 = ((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))^2*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))^2+(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(-m42+2*A2*(-m3+m1)+2*A1*(A2-m3+m2)+(-m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2))*(-m42+2*A2*(m3+m1)+2*A1*(A2+m3+m2)+(m3+m1+m2)^2-2*p1*p2*(ct1ct2+ch3h4*st1st2)) 

    C2::Float32 =-4*((p1*ct3ct1+p2*ct3ct2)+(p1*ch1h3*st3st1+p2*ch1h4*st3st2))*(m32-m42+2*A2*m1+2*A1*(A2+m2)+(m1+m2)^2-2*p1*p2*(ct1*ct2+ch3h4*st1st2))

    C4::Float32 = -8*(A1+A2+m1+m2+(p1*ct3ct1+p2*ct3ct2)+p1*ch1h3*st3st1+p2*ch1h4*st3st2)*(A1+A2+m1+m2-(p1*ct3ct1+p2*ct3ct2)-(p1*ch1h3*st3st1+p2*ch1h4*st3st2))

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


#=Testing
p1v = [1f0, 0.5f0, 1.8f0,]
p2v = [2f0, 0.2f0, 0.7f0]
p3v = [0f0 0f0; 0.3f0 0.3f0; 0.7f0 0.7f0]
p3v = [0f0;0.3f0; 0.7f0]
p3pv = copy(p3v)
p3v2 = [p3v p3pv]

Momentum3Value!(p3v2,p1v,p2v)
p3v2

Momentum3Value2!(p3v,p3pv,p1v,p2v)
p3v
p3pv

Momentum3Value3!(p3v,p3pv,p1v,p2v)
p3v
p3pv

p3v = Float32[550.49805, 0.6716571, 1.1593782]
p1v = Float32[1353.1063, 0.9887152, 1.2448756]
p2v = Float32[279.22748, 0.87663436, 1.7217331]

p3v = Float32[0.006277652, -0.43729234, 0.62806463]
p1v = Float32[0.008336283, -0.23659337, 0.636328]
p2v = Float32[0.0055497563, -0.08499956, 0.80733776]

p3v = Float32[7.147129, 0.7736238, 1.9076179]
p1v = Float32[9.2896595, 0.75850165, 1.802476]
p2v = Float32[5.56435, 0.6841588, 1.751785]

p3v = Float32[2399.8464, -0.9112277, 0.6679145]
p1v = Float32[4614.545, -0.737568, 0.6580323]
p2v = Float32[1833.5154, -0.98802865, 0.08539677]

p3v = Float32[2477.5803, 0.70000625, 1.2908807]
p2v = Float32[4426.429, 0.68573844, 1.2220751]
p1v = Float32[3568.4944, 0.15726769, 1.4209498]

p3v = Float32[550.49805, 0.6716571, 1.1593782]
p1v = Float32[1353.1063, 0.9887152, 1.2448756]
p2v = Float32[279.22748, 0.87663436, 1.7217331]

Sval = SValue(p3v,p1v,p2v)=#