function sCheck(sSmol::Float32,sBig::Float32)
    #returns true/false depending on if s is above minimum value
    s = sSmol+sBig
    return (s>=(mu1+mu2)^2)&&(s>=(mu3+mu4)^2) 

end

function tCheck(tSmol::Float32,tBig::Float32)
    #returns true/false depending on if t is above minimum value
    t = tSmol+tBig
    return (t<=(mu1-mu3)^2)&&(t<=(mu2-mu4)^2)

end

function uCheck(uSmol::Float32,uBig::Float32)
    #returns true/false depending on if t is above minimum value
    u = uSmol+uBig
    return (u<=(mu1-mu4)^2)&&(u<=(mu2-mu3)^2)

end

function stuCheck(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
    s = sSmol+sBig
    t = tSmol+tBig
    u = uSmol+uBig

    h = mu1^2+mu2^2+mu3^2+mu4^2
    a = (mu1*mu2-mu3*mu4)*(mu1^2+mu2^2-mu3^2-mu4^2)/h
    b = (mu1*mu3-mu2*mu4)*(mu1^2+mu3^2-mu2^2-mu4^2)/h
    c = (mu1*mu4-mu3*mu2)*(mu1^2+mu4^2-mu3^2-mu2^2)/h

    return (s*t*u>a*s+b*t+c*u)

end
