"""
    location(u,l,num,val)

Returns the index of the bin in which 'val' is contatined based on the 'num' of bins and their 'u' upper and 'l' lower bound.

# Examples
```julia-repl
julia> location(0f0,10f0,9,2f0)
2
```
"""
function location(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    return val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
end

"""
    locationp3(u,l,num,val)

Returns the index of the bin in which 'val' is contatined based on the 'num' of bins and their 'u' upper and 'l' lower bound including overflow and underflow possibilities. Index 1 reserved for underflow and index 2 for overflow, other bin indeicies are shifted by 2 i.e.  bin+2.

# Examples
```julia-repl
julia> locationp3(0f0,10f0,9,2f0)
4
```
"""
function locationp3(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    loc = val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
    return 1 <= loc <= num ? loc+2 : loc>num ? 2 : 1 # assignes 1 for under, 2 for over and loc+2 for in range
end
