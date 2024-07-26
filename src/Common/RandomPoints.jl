#=
Contains functions for generating random points
=#

"""
    RPointLogMomentum!(pu,pl,pv,num)

Edits the first element of `pv` with a random real-space momentum value between ``10^{pl}`` and ``10^{pu}``. This sample is chosen by first randomly picking a momentum bin in the range `1:num` and then uniformly sampling a momentum point in real-space (rather than log10 space) between pl and pu which are the momentum values at start and end of that bin. Sampling is done such there will be a constant number of points per momentum-space volume. As the momentum space between ``10^{pl}`` and ``10^{pu}`` it is a spherical shell hence the correct sampling is ``p = (U*(10^{pu})^3+(1-U)*(10^{pl})^3)^{1/3}`` with uniform ``U ∈ [0~~1]``.

Assumes ``f(x,p,μ)=f(x,\\vec{p})*(2πp^2)=const`` in bin, therefore momentum space volume element is ``\\mathrm{d}p`` and as such uniform sampling corresponds to ``U*10^{u}+(1-U)*10^{l}`` where ``U`` is a uniform random number between 0 and 1.

If instead ``f(x,\\vec{p})=const`` in bin, momentum space volume element is ``p^2 \\mathrm{d}p`` and uniform sampling corresponds to ``(10^pu)*\\sqrt[3]{U+(1-U)*10^{3pl-3pu}}`` where ``U`` is a uniform random number between 0 and 1.
"""
function RPointLogMomentum!(pv::Vector{Float64},pu::Float64,pl::Float64,num::Int64) 
    # Inputs a momentum vector and mometum bounds and mutates first of said vector
    bin = rand(1:num)
    l = (pl + (pu-pl)*(bin-1)/num)
    u = (pl + (pu-pl)*(bin)/num)

    U = rand(Float64)
    #pv[1] = (10^u)*cbrt(U+(1-U)*1f3^(l-u)) 
    pv[1] = U*10^(u)+(1-U)*10^(l)  # if instead want to sample space uniformly.

    return nothing
    
end

"""
    RPointSphereThetaPhi!()

Assigns the second (cos(theta)) and third (phi) elements of 'a' with a randomly, uniformly sampled values of spherical angles cos(theta) and phi (phi normalised by pi). 
"""
function RPointSphereCosThetaPhi!(a::Vector{Float64}) 
    # Inputs a 3 element vector [p, cos(theta), phi] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places)
    # phi points are normalised by pi

    u::Float64 = rand(Float64)
    v::Float64 = rand(Float64)

    a[2] = 2*v-1     # cos(theta) bound by [1,-1]
    a[3] = 2*u       # phi bound by [0,2) 

    return nothing
    
end