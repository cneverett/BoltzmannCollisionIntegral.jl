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

    n = abs((mEle^2*c^2)/(z2*ħ*q*B)) * p1 * (E2-p2*ct1*ct2)
    #println(n)

    y = p2 * st2 *st1 / (E2-p2*ct2*ct1) # y=x/n
    #println(x)

    # characteristic frequency
    #ω0 = abs((z2*q*B))/(E2*mEle)
    #println("critical photon momentum: "*string(ħ*ω0/(mEle*c^2)*E2^3))


    if n > 1e2 #&& 1-y < 0.01
        # approximation for J's to second order 
        e = 1-y^2
        K13 = besselk(1/3,n*e^(3/2)/3)
        K23 = besselk(2/3,n*e^(3/2)/3)
        J1 = ((sqrt(e))/(pi*sqrt(3)))*(K13 #=+(e/10)*(K13-2*n*e^(3/2)*K23)=#)
        J2 = (e/(pi*sqrt(3)))*(K23 #=+ (e/5)*(2*K23-(1/(e^(3/2)*n)+n*e^(3/2))*K13)=#)
    elseif n < 1e0
        # omega < omega0 therefore no synchrotron radiation
        J1 = 0e0
        J2 = 0e0
    else
        # exact J's
        J1 = besselj(n,n*y)
        J2 = 1/2 * (besselj(n-1,n*y) - besselj(n+1,n*y))
    end

    val = (abs(z2*m2^3/B))*(p1/E2)*((Jfactor1*J1)^2+(Jfactor2*J2)^2)
    #println(val)

    factor = (3*c^4*mEle^5)/(4*pi*ħ^3*μ0*q^3) # synchrotron emission rate divided by c*σT

    #println(factor)

    return val*factor
    
end

#SyncKernel([1e-14,0.6],[1e1,0.5],1,1,1e-4)
#n = 7.36e10
#e = 1-(0.0045)^2
#besselk(1/3,n*e^(3/2)/3)