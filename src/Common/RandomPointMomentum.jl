
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
    pv[1] = rand(Float32)*(pu-pl)+pl  # normal distribution

    return nothing
    
end

"""
    RPointLogMomentum!(pu,pl,pv)

Edits the first element of 'pv' with a random real-space momentum value by uniformly sampling Log10 space between 'pl' and 'pu' (which are in Log10 space). 
"""
function RPointLogMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32}) 
    # Inputs a momentum vector and mometum bounds and mutates first of said vector

    # generates random value uniform in log space
    pv[1] = 1f1^(rand(Float32)*(pu-pl)+pl)  # normal distribution

    return nothing
    
end


