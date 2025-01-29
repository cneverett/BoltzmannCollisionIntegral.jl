# =============== Values on Grid Boundaries ====================== #

"""
    bounds(up_bound,low_bound,num,spacing)

Returns a `num+1` long `Vector{Float}` of grid bounds. These grid bounds can spaced either by 
    - linear spacing: `spacing = "u"`
    - log10 spacing: `spacing = "l"`
    - binary (1/2^n) spacing: `spacing = "b"`

"""
function bounds(low_bound::T,up_bound::T,num::Int64,spacing::String) where T <: Union{Float32,Float64}

    if spacing == "u" # uniform spacing
        return [range(low_bound,up_bound,num+1);]
    elseif spacing == "l" # log spacing
        return 10 .^[range(low_bound,up_bound,num+1);]
    elseif spacing == "b" # binary (2^n) fractional spacing
        pow = [range(1,(num-1)/2); Inf]
        a = 1 .-(1/2) .^(pow)
        return [reverse(-a) ; a]
    else
        error("Spacing type not recognized")
    end
end

"""
    bounds_p(pl,pu,nump)

Returns a `nump+1` long `Vector{Float}` of p-space grid bounds NOT in Log10 space.

# Examples
```julia-repl
julia> bounds_p(-5e0,4e0,9)
10-element Vector{Float64}:
 1.0e-5
 1.0e-4
 1.0e-3
 0.01
 0.1
 1.0
 10.0
 100.0
 1000.0
 10000.0
```
"""
function bounds_p(pl::T,pu::T,nump::Int64) where T <: Union{Float32,Float64}
    # returns a vector{T} of p grid bounds NOT in Log10 space
    return 10 .^[range(pl,pu,nump+1);]
end

"""
    bounds_t(numt)

Returns a `numt+1` long `Vector{Float}` of theta-space grid bounds in terms of cos(theta).
Upper and lower bounds [tl tu] are defined as CONST in Init.jl as [-1 1], type returned is that of tl, tu.

# Examples
```julia-repl
julia> bounds_t(8)
9-element Vector{Float64}:
 -1.0
 -0.75
 -0.5
 -0.25
  0.0
  0.25
  0.5
  0.75
  1.0
```
"""
function bounds_t(numt::Int64) 
    return [range(tl,tu,numt+1);]
end

# ================================================================ #

# =============== Momentum and Angle Changes over grid cells ===== #
"""
    deltaVector(valr)

Inputs a `num+1` long `Vector{Float}` quantity values (domain bounds) and returns a `num` long `Vector{Float}` of differences (domain widths).

# Examples
```julia-repl
julia> deltaVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0])
3-element Vector{Float64}:
 9.0
 90.0
 900.0
```
"""
function deltaVector(valr::Vector{T}) where T <: Union{Float32,Float64}
    num = size(valr)[1]-1  # number of grid cells
    Δ = zeros(T,num)
    
    for ii in 1:num 
        Δ[ii] += abs(valr[ii+1] - valr[ii]) # abs accounts for Δμ = Δcosθ = cospi(t i) - cospi(t 1+1)
    end

    return Δ
end

# ================================================================ #

# ============== Mean Values in Grid Cells ======================= #
"""
    meanVector(valr)

Inputs a `num+1` long `Vector{Float}` of domain bounds and returns a `num` long `Vector{Float}` of mean value in domain range.

# Examples
```julia-repl
julia> meanVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0])
3-element Vector{Float64}:
 5.5
 55.0
 550.0
```
"""
function meanVector(valr::Vector{T}) where T <: Union{Float32,Float64}
    num = size(valr)[1]-1  # number of grid cells
    mean = zeros(T,num)
    
    for ii in 1:num 
        mean[ii] = (valr[ii+1] + valr[ii])/2
    end

    return mean
end

# ================================================================ #

# ========================= "Delta Energy" ======================= #
"""
    deltaEVector(pr,mu)

Inputs a `num+1` long `Vector{Float}` of p grid boundaries and the particle `mu` value (normalised mass) and returns a `num` long `Vector{Float}` of average energy values per grid cell.

# Examples
```julia-repl
julia> deltaEVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0], 1.0e0)
3-element Vector{Float64}:
 50.600693
 4951.15
 495001.16
```
"""
function deltaEVector(pr::Vector{T},mu::T) where T <: Union{Float32,Float64}

    num = size(pr)[1]-1  # number of grid cells
    E = zeros(T,num+1)
    ΔE = zeros(T,num)

    if mu == zero(T)
        for ii in 1:num 
            ΔE[ii] = pr[ii+1]^2-pr[ii]^2
            ΔE[ii] /= 2
        end
    else 
        for ii in 1:num+1 
            E[ii] = mu + pr[ii]^2/(sqrt(mu^2+pr[ii]^2)+mu)
        end 
        for ii in 1:num 
            ΔE[ii] = (pr[ii+1]-pr[ii])*mu
            ΔE[ii] += pr[ii+1]^3/(E[ii+1]+mu) - pr[ii]^3/(E[ii]+mu) 
            ΔE[ii] += mu^2*(asinh(pr[ii+1]/mu)-asinh(pr[ii]/mu))
            ΔE[ii] /= 2
        end
    end 

    return ΔE
end

# ================================================================ #

# ====================== "Delta Kinetic Energy" ================== #
"""
    deltaEkinVector(pr,mu)

Inputs a `num+1` long `Vector{Float}` of p grid boundaries and the particle `mu` value (normalised mass) and returns a `num` long `Vector{Float}` of average kinetic energy values per grid cell.

# Examples
```julia-repl
julia> deltaEkinVector([1.0e0, 10.0e0, 100.0e0, 1000.0e0], 1.0e0)
3-element Vector{Float64}:
     46.10069600605712
   4906.1506753523645
 494551.15128635924
```
"""
function deltaEkinVector(pr::Vector{T},mu::T) where T <: Union{Float32,Float64}
    # inputs a (num+1) vector{Float} of p grid boundaries and the particle mu value and return a (num) vector{Float} of average energy values per grid cell

    num = size(pr)[1]-1  # number of grid cells
    E = zeros(T,num+1)
    ΔEkin = zeros(T,num)

    if mu == zero(T) # same as ΔE
        for ii in 1:num 
            ΔEkin[ii] = pr[ii+1]^2-pr[ii]^2
            ΔEkin[ii] /= 2
        end
    else 
        for ii in 1:num+1 
            E[ii] = mu + pr[ii]^2/(sqrt(mu^2+pr[ii]^2)+mu)
        end 
        for ii in 1:num 
            # ΔEkin = ΔE - Δp*m
            ΔEkin[ii] = pr[ii+1]^3/(E[ii+1]+mu) - pr[ii]^3/(E[ii]+mu) 
            ΔEkin[ii] += mu^2*(asinh(pr[ii+1]/mu)-asinh(pr[ii]/mu))
            ΔEkin[ii] /= 2
        end
    end 

    return ΔEkin
end

# =============================================================== #