"""
    sCheck(sSmol,sBig,mu3,mu4)

Returns 'true' if 's' mandelstram generated from inital system state can generate a physical output state. 
"""
function sCheck(sSmol::Float64,sBig::Float64,mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)
    #returns true/false depending on if s is above minimum value
    # s >= (m1+m2)^2 && s >= (m3+m4)^2
    # sBig = (m1+m2)^2 in code
    return ((sSmol>(mu1+mu2)^2-sBig) && (sSmol>(mu3+mu4)^2-sBig))  # = gives T value of zero

end

function sCheck(sSmol::Float64,sBig::Float64,P::BinaryParameters)
    #returns true/false depending on if s is above minimum value
    mu1 = P.mu1
    mu2 = P.mu2
    mu3 = P.mu3
    mu4 = P.mu4
    # s >= (m1+m2)^2 && s >= (m3+m4)^2
    # sBig = (m1+m2)^2 in code
    return ((sSmol>(mu1+mu2)^2-sBig) && (sSmol>(mu3+mu4)^2-sBig))  # = gives T value of zero

end

function tCheck(tSmol::Float64,tBig::Float64,mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)
    #returns true/false depending on if t is above minimum value
    # t <= (m1-m3)^2 && t <= (m2-m4)^2
    # tBig = (m3-m1)^2 in code
    return ((tSmol<(mu1-mu3)^2-tBig) && (tSmol<(mu2-mu4)^2-tBig))

end

function uCheck(uSmol::Float64,uBig::Float64,mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)
    #returns true/false depending on if t is above minimum value
    # u <= (m1-m4)^2 && u <= (m2-m3)^2
    # uBig = (m2-m3)^2 in code
    return ((uSmol<(mu1-mu4)^2-uBig) && (uSmol<(mu2-mu3)^2-uBig))

end

function stuCheck(sSmol::Float64,sBig::Float64,tSmol::Float64,tBig::Float64,uSmol::Float64,uBig::Float64,mu1::Float64,mu2::Float64,mu3::Float64,mu4::Float64)
    s = sSmol+sBig
    t = tSmol+tBig
    u = uSmol+uBig

    h = mu1^2+mu2^2+mu3^2+mu4^2
    a = (mu1*mu2-mu3*mu4)*(mu1^2+mu2^2-mu3^2-mu4^2)/h
    b = (mu1*mu3-mu2*mu4)*(mu1^2+mu3^2-mu2^2-mu4^2)/h
    c = (mu1*mu4-mu3*mu2)*(mu1^2+mu4^2-mu3^2-mu2^2)/h

    check2::Bool = s*t*u >= a*s+b*t+c*u

    if check2 == false
        println(s+t+u)
    end

    return check2

end
