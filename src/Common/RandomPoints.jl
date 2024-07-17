#=
Contains functions for generating random points
=#

"""
    RPointLogMomentum!(pu,pl,pv,num)

Edits the first element of 'pv' with a random real-space momentum value between 10^pl and 10^pu. This sample is chosen by first randomly picking a momentum bin in the range 1:num and then uniformly sampling a momentum point in real-space (rather than log10 space) between l and u which are the momentum values at start and end of that bin. Sampling is done such there will be a constant number of points per momentum-space volume. As the momentum space between 10^pl and 10^pu it is a spherical shell hence the correct sampling is p = (U*(10^pu)^3+(1-U)*(10^pl)^3)^1/3 with uniform U ∈ [0 1].

Assumes ``f(x,p,μ)=f(x,\\vec{p})*(2pi p^2)=const`` in bin, therefore momentum space volume element is ``dp`` and as such uniform sampling corresponds to ``U*10^(u)+(1-U)*10^(l)`` where U is a uniform random number between 0 and 1.

If instead ``f(x,\\vec{p})=const`` in bin, momentum space volume element is ``p^2 dp`` and uniform sampling corresponds to ``(10^u)*cbrt(U+(1-U)*10^(3l-3u))`` where U is a uniform random number between 0 and 1.
"""
function RPointLogMomentum!(pu::Float32,pl::Float32,pv::Vector{Float32},num::Int64) 
    # Inputs a momentum vector and mometum bounds and mutates first of said vector
    bin = rand(1:num)
    l = (pl + (pu-pl)*(bin-1)/num)
    u = (pl + (pu-pl)*(bin)/num)

    U = rand(Float32)
    #pv[1] = (1f1^u)*cbrt(U+(1f0-U)*1f3^(l-u)) 
    pv[1] = U*1f1^(u)+(1f0-U)*10^(l)  # if instead want to sample space uniformly.

    return nothing
    
end

"""
    RPointSphereThetaPhi!()

Assigns the second (cos(theta)) and third (phi) elements of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta) and phi (phi normalised by pi). 
"""
function RPointSphereCosThetaPhi!(a::Vector{Float32}) 
    # Inputs a 3 element vector [p, cos(theta), phi] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places)
    # phi points are normalised by pi

    u::Float32 = rand(Float32)
    v::Float32 = rand(Float32)

    a[2] = 2*v-1     # cos(theta) bound by [1,-1]
    a[3] = 2*u       # phi bound by [0,2) 

    return nothing
    
end