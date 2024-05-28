#=

This module defines the differential/total cross section function and its normalisation for specific interactions, it therefore requires the particle names to be defined as global constants for use in determing the correct cross sections functions to compile.

=#

#= if (name1 == name2 == name3 == name4) # all particles identical (1S 1T)
    export dsigmadt, dsigmadtn, sigma, sigman
elseif ((name1 == name2) && (name3 == name4))   # input particles are identical and final particles are identical (2T 2S)
    error("not yet implimented")
elseif ((name1 == name3) && (name2 == name4))   # input states and final state are identical (1T 2S)
    error("not yet implimented")
elseif (name1 == name2) && (name3 != name4) ||  (name1 != name2) && (name3 == name4) # either states 12 (or 34) are identical, but the other states 34 (or 12)  are not ()
    error("not yet implimented")
elseif (name1 != name2 != name3 != name4) # not particles identical
    error("not yet implimented")
else
    error("particle states not considered")
end =#

# Dependancies  
include("MyPhysicalConstants.jl")

#======================= Hard Sphere Collisions ===============================#

if (name1 == "Sph" && name2 == "Sph" && name3 == "Sph" && name4 == "Sph")
    # Hard sphere collisions
    function dsigmadt(s::Float32,t::Float32)
        
        1f0/(s-4*muSph^2)

    end

    function dsigmadt(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
        
        1f0/(sSmol#=+(sBig-4*muSph^2)=#) # sBig = (m3+m4)^2=4musph^2

    end

    function sigma(s::Float32)
        
        1f0/2f0 # factor of 2 accounts for identical final states

    end

    function sigma(sSmol::Float32,sBig::Float32)
        
        1f0/2f0 # factor of 2 accounts for identical final states

    end

    const dsigmadtn = Float32(pi)*(2f0*RSph)^2
    const sigman = Float32(pi)*(2f0*RSph)^2

    return nothing

end

#=======================================================================#

#============== Electron Positron Annihilation to Two Photons ==========#

if (name1 == "Ele" && name2 == "Pos" && name3 == "Pho" && name4 == "Pho")
    # Hard sphere collisions
    function dsigmadt(s::Float32,t::Float32) # Berestetskii (88.4)
        
        -(1/(s(s-4)))*((1/(t-1)+1/(1-s-t))^2+(1/(t-1)+1/(1-s-t))-(1/4)*((t-1)/(1-s-t)+(1-s-t)/(t-1)))

    end

    function dsigmadt(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
        
        s = sSmol+sBig # sBig = (m1+m2)^2 = 4 (normalised units) -> s = sSmol + 4
        #t = tSmol+tBig # tBig = (m3-m1)^2 = 1 (normalised units) -> t = tSmol + 1
        #u = uSmol+uBig # uBig = (m2-m3)^2 = 1 (normalised units) -> u = uSmol + 1
        -(1/((s)*(sSmol)))*((1/(tSmol)+1/(uSmol))^2+(1/(tSmol)+1/(uSmol))-(1/4)*((tSmol)/(uSmol)+(uSmol)/(tSmol)))

    end

    function sigma(s::Float32) # Berestetskii (88.6)
        
        (1/(4*s^2*(s-4)))*((s^2+4*s-8)*log((sqrt(s)+sqrt(s-4))/(sqrt(s)-sqrt(s-4)))-(s+4)*sqrt(s*(s-4)))

    end

    function sigma(sSmol::Float32,sBig::Float32)

        s = sSmol+sBig
        (1/(4*(sSmol)*s^2))*((sSmol^2+12*sSmol+24)*log((s+sSmol+2*sqrt(sSmol*s))/(sBig))-(sSmol+8)*sqrt((s)*(sSmol)))

    end

    const dsigmadtn = 3*σT;
    const sigman = 3*σT;

    return nothing

end

# ==================================================================== # 


#======== Electron Positron Pair Production from Two Photons ==========#

if (name1 == "Pho" && name2 == "Pho" && name3 == "Ele" && name4 == "Pos")
    # Hard sphere collisions
    function dsigmadt(s::Float32,t::Float32)
        
        -(1/(s^2))*((1/(t-1)+1/(1-s-t))^2+(1/(t-1)+1/(1-s-t))-(1/4)*((t-1)/(1-s-t)+(1-s-t)/(t-1)))

    end

    function dsigmadt(sSmol::Float32,sBig::Float32,tSmol::Float32,tBig::Float32,uSmol::Float32,uBig::Float32)
        
        s = sSmol # sBig = (m1+m2)^2 = 0 (normalised units) -> s = sSmol
        #t = tSmol+tBig # tBig = (m3-m1)^2 = 1 (normalised units) -> t = tSmol + 1
        #u = uSmol+uBig # uBig = (m2-m3)^2 = 1 (normalised units) -> u = uSmol + 1
        -(1/(sSmol^2))*((1/(tSmol)+1/(uSmol))^2+(1/(tSmol)+1/(uSmol))-(1/4)*((tSmol)/(uSmol)+(uSmol)/(tSmol)))

    end

    function sigma(s::Float32) 
        
        (1/(2*s^3))*((s^2+4*s-8)*log((sqrt(s)+sqrt(s-4))/(sqrt(s)-sqrt(s-4)))-(s+4)*sqrt(s*(s-4)))

    end

    function sigma(sSmol::Float32,sBig::Float32)

        # For photon-photon annihilation, sBig=0 and sSmol=s but still want to avoid float issues with s-4
        s = sSmol+sBig
        (1/(2*s^3))*((s^2+4*s-8)*log((2*s-4+2*sqrt(s*(s-4)))/(4))-(s+4)*sqrt(s*(s-4)))

    end

    const dsigmadtn = 3*σT;
    const sigman = 3*σT;

    return nothing

end

# ==================================================================== # 