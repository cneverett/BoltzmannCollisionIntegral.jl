#=

This module generates vectors of useful grid quantities from grid upper, lower bounds and number of grid cells.
The "useful quantities" are then used in calculating moments of the distribution function.

    0-th Moment:
    ∫f(ij)sin(θ)dpdθ = ∑(ij) fij Δpi Δμj
    where:
        p(i+1/2) are the momentum values on grid boundries
        θ(j+1/2) are the theta values on grid boundries 
        Δpi = p(i+1/2) - p(i-1/2)
        Δμj = cos(θ(j-1/2)) - cos(θ(j+1/2)) = -[cos(θ(j+1/2)) - cos(θ(j-1/2))]
    all are positive definite

    1-st Moment:
    ∫p^μ f(ij)sin(θ) dpdθ
        Momentum part: 
        ∫p f(ij) sin(θ) dp dθ = ∑(ij) f(ij) Δμj Δpi <p>i
        where:
            <p>i = 1/2(p(i+1/2)-p(i-1/2)) is the mean momnetum
=#

# =============== Values on Grid Boundaries ====================== #
"""
    prange(pl,pu,nump)

Returns a `nump+1` long `Vector{Float}` of p-space grid bounds NOT in Log10 space.

# Examples
```julia-repl
julia> prange(-5f0,4f0,9)
10-element Vector{Float32}:
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
function prange(pl::T,pu::T,nump::Int64) where T <: Union{Float32,Float64}
    # returns a vector{T} of p grid bounds NOT in Log10 space
    return 10 .^[range(pl,pu,nump+1);]
end

"""
    trange(numt)

Returns a `numt+1` long `Vector{Float}` of theta-space grid bounds in terms of cos(theta).
Upper and lower bounds [tl tu] are defined as CONST in Init.jl as [-1 1], type returned is that of tl, tu.

# Examples
```julia-repl
julia> trange(8)
9-element Vector{Float32}:
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
function trange(numt::Int64) 
    return [range(tl,tu,numt+1);]
end

# ================================================================ #

# =============== Momentum and Angle Changes over grid cells ===== #

"""
    deltaVector(valr)

Inputs a `num+1` long `Vector{Float}` quantitiy values (domain bounds) and returns a `num` long `Vector{Float}` of differeces (domain widths).

# Examples
```julia-repl
julia> deltaVector([1.0f0, 10.0f0, 100.0f0, 1000.0f0])
3-element Vector{Float32}:
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
julia> meanVector([1.0f0, 10.0f0, 100.0f0, 1000.0f0])
3-element Vector{Float32}:
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

Inputs a `num+1` long `Vector{Float}` of p grid boundries and the particle `mu` value (normalised mass) and returns a `num` long `Vector{Float}` of average energy values per grid cell.

# Examples
```julia-repl
julia> deltaEVector([1.0f0, 10.0f0, 100.0f0, 1000.0f0], 1.0f0)
3-element Vector{Float32}:
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
