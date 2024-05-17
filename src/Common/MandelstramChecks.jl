function sCheck(sSmol::Float32,sBig::Float32)
    #returns true/false depending on if s is above minimum value
    # s >= (m1+m2)^2 && s >= (m3+m4)^2
    # sBig = (m1+m2)^2 in code
    return (sSmol>0f0)&&(sSmol>(mu3+mu4)^2-sBig)  # = gives T value of zero

end

function tCheck(tSmol::Float32,tBig::Float32)
    #returns true/false depending on if t is above minimum value
    # t <= (m1-m3)^2 && t <= (m2-m4)^2
    # tBig = (m3-m1)^2 in code
    t = tSmol+tBig
    return (tSmol<0f0)&&(tSmol<(mu2-mu4)^2-tBig)

end

function uCheck(uSmol::Float32,uBig::Float32)
    #returns true/false depending on if t is above minimum value
    # u <= (m1-m4)^2 && u <= (m2-m3)^2
    # uBig = (m2-m3)^2 in code
    u = uSmol+uBig
    return (uSmol<(mu1-mu4)^2-uBig)&&(uSmol<0f0)

end

function stuCheck(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
    s = sSmol+sBig
    t = tSmol+tBig
    u = uSmol+uBig

    h = mu1^2+mu2^2+mu3^2+mu4^2
    a = (mu1*mu2-mu3*mu4)*(mu1^2+mu2^2-mu3^2-mu4^2)/h
    b = (mu1*mu3-mu2*mu4)*(mu1^2+mu3^2-mu2^2-mu4^2)/h
    c = (mu1*mu4-mu3*mu2)*(mu1^2+mu4^2-mu3^2-mu2^2)/h

    check2 = s*t*(h-s-t) >= a*s+b*t+c*(h-s-t)

    return check2

end
