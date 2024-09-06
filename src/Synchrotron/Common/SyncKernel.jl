function SyncKernel(p1v,p2v,m2,z2,B)

    # p1 is Photon
    # p2 is Charged Particle

    p1 = p1v[1]
    p2 = p2v[1]
    ct1 = p1v[2]
    ct2 = p2v[2]
    st1 = sqrt(1-ct1^2)
    st2 = sqrt(1-ct2^2)
    
    E2 = sqrt(p2^2 + m2^2)

    Jfactor1 = (E2*ct1-p2*ct1*ct2)/(st1)
    Jfactor2 = p2*st2

    n = (mEle^2*c^2)/(z2*ħ*q*B) * p1 * (E2-p2*ct1*ct2)
    println(n)

    x = p2 * st2 *st1 / (E2-p2*ct2*ct1)
    println(x)
    e = 1-x^2

    # characteristic frequency
    ω0 = (z2*q*B)/(E2*mEle)
    println("critical photon momentum: "*string(ħ*ω0/(mEle*c^2)*E2^3))


    if n > 1e4 && 1-x < 0.001
        # approximation for J's to second order 
        K13 = besselk(1/3,n*e^(3/2)/3)
        K23 = besselk(2/3,n*e^(3/2)/3)
        J1 = ((sqrt(e))/(pi*sqrt(3)))*(K13 +(e/10)*(K13-2*n*e^(3/2)*K23))
        J2 = (e/(pi*sqrt(3)))*(K23 + (e/5)*(2*K23-(1/(e^(3/2)*n)+n*e^(3/2))*K13))
    else
        # exact J's
        J1 = besselj(n,n*x)
        J2 = 1/2 * (besselj(n-1,x) - besselj(n+1,x))
    end

    val = (z2/B)*(p1^3/E2)*((Jfactor1*J1)^2+(Jfactor2*J2)^2)

    return val
    
end
