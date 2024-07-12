
#= 
This script contains allocating and non-allocating functions for generating uniformly distributed, random points on the surface of a unit sphere and half sphere 
=#

"""
    RPointMomentum!(pu,pl,pv)

Edits the first element of 'pv' with a random, uniformly sampled, momentum value between 'pl' and 'pu'
"""
function RPointMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32}) 
    # Inputs a momentum vector and mometum bounds and mutates firstof said vector

    # generates random value uniform in real space
    pv[1] = rand(Float32)*(pu-pl)+pl

    return nothing
    
end

"""
    RPointLogMomentum!(pu,pl,pv)

Edits the first element of 'pv' with a random real-space momentum value by uniformly sampling Log10 space between 'pl' and 'pu' (which are in Log10 space). 
"""
function RPointLogMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32}) 
    # Inputs a momentum vector and mometum bounds and mutates first of said vector

    # generates random value uniform in log space
    pv[1] = 1f1^(rand(Float32)*(pu-pl)+pl)  

    return nothing
    
end

"""
    RPointLogMomentum!(pu,pl,pv,num)

Edits the first element of 'pv' with a random real-space momentum value between 10^pl and 10^pu. This sample is chosen by first randomly picking a momentum bin in the range 1:num and then uniformly sampling a momentum point in real-space (rather than log10 space) between l and u which are the momentum values at start and end of that bin. Sampling is done such there will be a constant number of points per momentum-space volume. As the momentum space between 10^pl and 10^pu it is a spherical shell hence the correct sampling is p = (U*(10^pu)^3+(1-U)*(10^pl)^3)^1/3 with uniform U ∈ [0 1].  
"""
function RPointLogMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32},num::Int64) 
    # Inputs a momentum vector and mometum bounds and mutates first of said vector
    bin = rand(1:num)
    l = (pl + (pu-pl)*(bin-1)/num)
    u = (pl + (pu-pl)*(bin)/num)

    U = rand(Float32)
    #pv[1] = cbrt(U*1f1^(3*u)+(1f0-U)*1f1^(3*l)) 
    pv[1] = U*10^(u)+(1f0-U)*10^(l) 

    return nothing
    
end