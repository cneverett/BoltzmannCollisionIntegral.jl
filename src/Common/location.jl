function location(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    return val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
end

function locationp3(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array. Bins MUST be uniform
    loc = val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
    return 1 <= loc <= num ? loc+2 : loc>num ? 2 : 1 # assignes 1 for under, 2 for over and loc+2 for in range
end
