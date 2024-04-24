    
function Momentum3Value!(p3v::Array{Float32},p1v::Vector{Float32},p2v::Vector{Float32},m1::Float32,m2::Float32,m3::Float32,m4::Float32)

    # pv should be [p,t,h]

    p1::Float32 = p1v[1]
    p2::Float32 = p2v[1]

    ct3::Float32 = cospi(p3v[2,1]) # sinpi and cospi slightly slower than sin(pi*) but more accurate apparently
    ct1::Float32 = cospi(p1v[2])
    ct2::Float32 = cospi(p2v[2]) 

    st3::Float32 = sinpi(p3v[2,1])
    st1::Float32 = sinpi(p1v[2])
    st2::Float32 = sinpi(p2v[2])

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

    # val = (C2-C3)/C4
    # valp =(C2+C3)/C4

    p3v[1,1] = 0f0 # reset p3v values
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
            if val >= 0f0 # parallel
                p3v[1,1] = val 
            else #if valp < 0f0  # anti-parallel 
                p3v[1,1] = -val 
                p3v[2,1] = mod(1f0-p3v[2,1],1f0)     # theta bound by [0,1]
                p3v[3,1] = mod(p3v[3,1]+1f0,2f0)     # phi bound by [0,2) 
            #else
                #error("p3 state not accounted for"*string(val))
            end
        else
            p3v[1,1] = 0f0
        end

        if (( valp != val) && (m1+m2-m3+p12/(sqm1p1+m1)+p22/(sqm2p2+m2)-valp^2/(sqrt(m32+valp^2)+m3))>0) # avoids counting same state twice
            #assign primed for aligned case as valp +ve
            if valp >= 0f0 # parallel
                p3v[1,2] = valp 
            else #if valp < 0f0  # anti-parallel 
                p3v[1,2] = -valp 
                p3v[2,2] = mod(1f0-p3v[2,2],1f0)     # theta bound by [0,1]
                p3v[3,2] = mod(p3v[3,2]+1f0,2f0)     # phi bound by [0,2) 
            #else
            # error("p3p state not accounted for:"*string(valp))
            end
        else
            p3v[1,2] = 0f0
        end

    else
        p3v[1,1] = 0f0
        p3v[1,2] = 0f0
    end

    return nothing

end



