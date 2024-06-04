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

Returns the index of the bin in which 'val' is contatined based on the 'num' of bins and their 'u' upper and 'l' lower bound including overflow and underflow possibilities. Overflow are assigned to num+1 while underflow are assigned to lowest bin i.e. 1.

# Examples
```julia-repl
julia> locationp3(10f0,1f0,9,2f0)
2
julia> locationp3(10f0,1f0,9,11f0) # overflow
10
julia> locationp3(10f0,1f0,9,0.5f0) # underflow
1
```
"""
function locationp3(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    loc = val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
    return 1 <= loc <= num ? loc : loc>num ? num+1 : 1 # assignes 1 for under, num+1 for over and loc for in range
end
