
#= 
This script contains allocating and non-allocating functions for generating uniformly distributed, random points on the surface of a unit sphere and half sphere 
=#


function RPointMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32}) 
    # Inputs a momentum vector and mometum bounds and mutates firstof said vector

    # generates random value uniform in real space
    pv[1] = rand(Float32)*(pu-pl)+pl  # normal distribution

    return nothing
    
end

function RPointLogMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32}) 
    # Inputs a momentum vector and mometum bounds and mutates first of said vector

    # generates random value uniform in log space
    pv[1] = 1f1^(rand(Float32)*(pu-pl)+pl)  # normal distribution

    return nothing
    
end


