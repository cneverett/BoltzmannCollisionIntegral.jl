#= 
This module contains allocating and non-allocating functions for generating uniformly distributed, random points on the surface of a unit sphere and half sphere 
=#

"""
    RPointSphere!(v)

Assigns 'v' with a random (x,y,z) point on the surface of a unit sphere.
"""
function RPointSphere!(v::Vector{Float32}) 
    # Inputs a vector and mutates said vector

    #nv = 0f0
    #while nv < 1f-2
    v[1] = randn(Float32)  # normal distribution
    v[2] = randn(Float32)
    v[3] = randn(Float32)
    nv::Float32 = sqrt(v[1]^2+v[2]^2+v[3]^2)
    #end

    v .= v ./ nv  

    return nothing
    
end

"""
    RPointSphere()

Returns a three element vecotr 'v' with a random (x,y,z) point on the surface of a unit sphere.
"""
function RPointSphere() 
   
    # Returns a vector
    v = Vector{Float32}(undef,3)
    #nv = 0f0
    #while nv < 1f-4
    v[1] = randn(Float32)  # normal distribution
    v[2] = randn(Float32)
    v[3] = randn(Float32)
    nv::Float32 = sqrt(v[1]^2+v[2]^2+v[3]^2)
    #end

    v .= v ./ nv  

    return v
    
end

"""
    RPointSphereThetaPhi!()

Assigns the second (theta) and third (phi) elements of 'a' with a randomly, uniformly sampled values of spherical angles theta and phi normalised by pi. 
"""
function RPointSphereThetaPhi!(a::Vector{Float32}) 
    # Inputs a 3 element vector [p, theta, phi] and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places) and results normalised by pi

    u::Float32 = rand(Float32)
    v::Float32 = rand(Float32)

    a[2] = acos(2*v-1) / pi     # theta bound by [0,1]
    a[3] = 2*u                  # phi bound by [0,2) 

    return nothing
    
end

"""
    RPointSphereThetaPhi!()

Assigns the second (theta) and third (phi) rows of array 'a' (3x2) with a randomly, uniformly sampled values of spherical angles theta and phi normalised by pi. 
"""
function R2PointSphereThetaPhi!(a::Array{Float32}) 
    # Inputs a 6 element vector ([p, theta, phi],[p', theta', phi']) and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places) and results normalised by pi
    # primed points are antiparallel to unprimed therefore theta'=mod(1-theta,1) & phi'=mod(phi+1,2)

    u::Float32 = rand(Float32)
    v::Float32 = rand(Float32)

    a[2,1]= acos(2*v-1) / pi     # theta bound by [0,1]
    a[3,1] = 2*u                  # phi bound by [0,2) 
    # primed angles
    a[2,2] = a[2,1]
    a[3,2] = a[3,1]

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

"""
    RPointSphereThetaPhi!()

Assigns the second (cos(theta)) and third (phi) rows of array 'a' (3x2) with a randomly, uniformly sampled values of spherical angles cos(theta) and phi (phi normalised by pi). 
"""
function R2PointSphereCosThetaPhi!(a::Array{Float32}) 
    # Inputs a 6 element vector ([p, theta, phi],[p', theta', phi']) and mutates said vector with new random values using form given in https://mathworld.wolfram.com/SpherePointPicking.html (with theta and phi changed places)
    # phi values are normalised by pi

    u::Float32 = rand(Float32)
    v::Float32 = rand(Float32)

    a[2,1] = 2*v-1    # cos(theta) bound by [1,-1]
    a[3,1] = 2*u     # phi bound by [0,2) 
    # primed angles
    a[2,2] = a[2,1]
    a[3,2] = a[3,1]

    return nothing
    
end

function RPointHalfSphere!(v::Vector{Float32}) 
    # Inputs a vector and mutates said vector, where the y value is always positive

    nv::Float32 = 0f0
    #while nv < 1f-4
    v[1] = randn()  # normal distribution
    v[2] = randn()
    v[3] = randn()
    nv = sqrt(v[1]^2+v[2]^2+v[3]^2)
    #end

    v .= v ./ nv;
    v[2] *= sign(v[2]) # makes y coordinate positive 

    return nothing
    
end

function RPointHalfSphere() 
    # returns a vector (allocating)
    
    v = Vector{Float32}(undef,3)
    #nv = 0f0
    #while nv < 1f-4
    v[1] = randn()  # normal distribution
    v[2] = randn()
    v[3] = randn()
    nv = sqrt(v[1]^2+v[2]^2+v[3]^2)
    #end

    v .= v ./ nv;
    v[2] *= sign(v[2]) # makes y coordinate positive 

    return v
    
end
