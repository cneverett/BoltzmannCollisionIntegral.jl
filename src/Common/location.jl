function location(u::Float32,l::Float32,num::Int64,val::Float32)
    # function for generating poisition in array
    return val != l ? ceil(Int64,Float32(num)*(val-l)/(u-l)) : Int64(1) 
end
