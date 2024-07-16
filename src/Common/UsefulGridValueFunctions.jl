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

    function prange(pl::Float32,pu::Float32,nump::Int64)
        # returns a vector{Float32} of p grid bounds NOT in Log10 space
        return 10 .^[range(pl,pu,nump+1);]
    end

    function trange(numt::Int64)
        # returns a vector{Float32} of t grid bounds in terms of cospi(t)
        # tl and tu defined as CONST in Init.jl 
        return [range(tl,tu,numt+1);]
    end

# ================================================================ #

# =============== Momentum and Angle Changes over grid cells ===== #

    function deltaVector(valr::Vector{Float32})
        # inputs a (num+1) vector{Float32} quantitiy values and returns a (num) vector{Float32} of differeces 

        num = size(valr)[1]-1  # number of grid cells
        Δ = zeros(Float32,num)
        
        for ii in 1:num 
            Δ[ii] += abs(valr[ii+1] - valr[ii]) # abs accounts for Δμ = Δcosθ = cospi(t i) - cospi(t 1+1)
        end

        return Δ
    end

# ================================================================ #

# ============== Mean Values in Grid Cells ======================= #

    function meanVector(valr::Vector{Float32})
        # inputs a (num+1) vector{Float32} quantitiy values and returns a (num) vector{Float32} of averages
        num = size(valr)[1]-1  # number of grid cells
        mean = zeros(Float32,num)
        
        for ii in 1:num 
            mean[ii] = (valr[ii+1] + valr[ii])/2
        end

        return mean
    end

# ================================================================ #

# ========================= "Delta Energy" ======================= #

    function deltaEVector(pr::Vector{Float32},mu::Float32)
        # inputs a (num+1) vector{Float32} of p grid boundries and the particle mu value and return a (num) vector{Float32} of average energy values per grid cell
        num = size(pr)[1]-1  # number of grid cells
        E = zeros(Float32,num+1)
        ΔE = zeros(Float32,num)

        for ii in 1:num+1 
            E[ii] = mu + pr[ii]^2/(sqrt(mu^2+pr[ii]^2)+mu)
        end 

        for ii in 1:num 
            ΔE[ii] = (pr[ii+1]-pr[ii])*mu
            ΔE[ii] += pr[ii+1]^3/(E[ii+1]+mu) - pr[ii]^3/(E[ii]+mu) 
            ΔE[ii] += mu^2*(asinh(pr[ii+1]/mu)-asinh(pr[ii]/mu))
            ΔE[ii] /= 2f0
         end 

        return ΔE
    end

# ================================================================ #
