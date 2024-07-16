"""
    location(u,l,num,val)

Returns the index of the bin in which 'val' is contatined based on the 'num' of bins and their 'u' upper and 'l' lower bound.

# Examples
```julia-repl
julia> location(10f0,0f0,9,2f0)
2
```
"""
function location(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    return val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
end

"""
    location_t(numt,val)

Returns the index of the bin in which the costheta 'val' is contatined based on the 'numt' of bins. Bounds [tl tu] are defined as CONST in Init.jl

# Examples
```julia-repl
julia> location_t(8,0.5f0)
6
```
"""
function location_t(numt::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    return val != tl ? ceil(Int64,Float32(numt)*(val-tl)/(tu-tl)) : Int64(1) 
end


"""
    location_p3(u,l,num,val)

Returns the index of the bin in which 'val' is contatined based on the 'num' of bins and their 'u' upper and 'l' lower bound including overflow and underflow possibilities. Overflow are assigned to num+1 while underflow are assigned to lowest bin i.e. 1.

# Examples
```julia-repl
julia> location_p3(10f0,1f0,9,2f0)
2
julia> location_p3(10f0,1f0,9,11f0) # overflow
10
julia> location_p3(10f0,1f0,9,0.5f0) # underflow
1
```
"""
function location_p3(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    logp = log10(val)
    loc = logp != l ? ceil(Int64,Float32(num)*(logp-l)/(u-l)) : Int64(1) 
    return 1 <= loc <= num ? loc : loc>num ? num+1 : 1 # assignes 1 for under, num+1 for over and loc for in range
end

"""
    vectorLocation(u1,l1,u2,l2,num1,num2,vector)

Returns a tuple of bin location for (log10momentum,cos(theta)) based on an input 'vector' and bounds 'u,l' of their domains and the 'num' of uniformly spaced bins.
costheta bounds [tl tu] are defined as CONST in Init.jl
"""
function vectorLocation(pu::Float32,pl::Float32,nump::Int64,numt::Int64,vector::Vector{Float32})

    logp = log10(vector[1])
    ctheta = vector[2]

    logploc = (logp != pl ? ceil(Int64,Float32(nump)*(logp-pl)/(pu-pl)) : Int64(1))
    cthetaloc = (ctheta != tu ? ceil(Int64,Float32(numt)*(ctheta-tl)/(tu-tl)) : Int64(1))
    
    return (logploc,cthetaloc)
    
end

"""
    vectorLocation(u1,l1,u2,l2,num1,num2,vector)

Returns a tuple of bin location for (log10momentum,cos(theta)) based on an input 'vector' and bounds 'u,l' of their domains and the 'num' of uniformly spaced bins, including overflow and underflow possibilities. Overflow are assigned to num1+1 while underflow are assigned to lowest bin i.e. 1. 
"""
function vectorLocationp3(u1::Float32,l1::Float32,u2::Float32,l2::Float32,num1::Int64,num2::Int64,vector::Vector{Float32})

    logp = log10(vector[1])
    ctheta = vector[2]

    logploc = (logp != l1 ? ceil(Int64,Float32(num1)*(logp-l1)/(u1-l1)) : Int64(1))
    cthetaloc = (ctheta != l2 ? ceil(Int64,Float32(num2)*(ctheta-l2)/(u2-l2)) : Int64(1))

    logploc = (1 <= logploc <= num1 ? logploc : logploc>num1 ? num1+1 : 1)
    
    return (logploc,cthetaloc)
    
end

#= # testing
p1v = [ 1.258741f5
0.93817377f0
 1.2125252f0]
 a=p1v[1]
@time location(10f0,0f0,9,p1v)

@time vectorLocation(10f0,1f0,1f0,-1f0,10,20,p1v)

@time vectorLocationp3(5f0,-4f0,1f0,-1f0,10,20,p1v)

locationp3(10f0,1f0,10,log10(a)) =#