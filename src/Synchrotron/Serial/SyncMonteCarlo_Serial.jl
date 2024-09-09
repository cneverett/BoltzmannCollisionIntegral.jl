"""
    SyncMonteCarloAxi_Serial!(SAtotal,SAtally,pMax,tMinMax,Parameters,numSiter)

# Arguments
- `SAtotal::Array{Float64,4}` : Array of stored integration totals for S matrix for 2+B->2+B+1 interaction
- `SAtally::Array{UInt32,4}` : Array of stored integration tallies for S matrix for 2+B->2+B+1 interaction
- `pMax::Array{Float64,3}` : Array of maximum momentum values for species 2
- `tMinMax::Array{Float64,3}` : Array of minimum and maximum theta values for species 2
- `Parameters::Tuple{Float64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numSiter::Int64` : Number of S iterations

# Output:
- Argument arrays SAtotal,SAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Breif
- Random Sample points in each of these domains
    - RandomPointSphere for theta (for species 1,2)
    - RandomPointMomentum for p ( species 1,2)
- Take random points (p1,p2,t1,t2) and calculate Synchrotron emissivity
- Find position in S arrays and allocated tallies and totals accordingly.
"""
function SyncMonteCarloAxi_Serial!(SAtotal::Array{Float64,4},SAtally::Array{UInt32,4},#=pMax::Array{Float64,5},tMinMax::Array{Float64,6}=#Parameters::Tuple{Float64,Float64,Float64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64},numTiter::Int64,numSiter::Int64)

    # Set Parameters
    (mu2,z2,BMag,p1l,p1u,nump1,p2l,p2u,nump2,numt1,numt2) = Parameters

    # allocate arrays
    p1v::Vector{Float64} = zeros(Float64,2)
    p2v::Vector{Float64} = zeros(Float64,2)
    Sval::Float64 = 0e0

    for _ in 1:numTiter

        # generate p2v (emitting particle)
        RPointSphereCosTheta!(p2v)
        RPointLogMomentum!(p2v,p2l,p2u,nump2)

        for _ in 1:numSiter

            # generate p1v (photon)
            RPointSphereCosTheta!(p1v)
            RPointLogMomentum!(p1v,p1l,p1u,nump1)

            # calculate S value
            Sval = SyncKernel(p1v,p2v,mu2,z2,BMag)
            # find S array location 
            (p1loc,t1loc) = vectorLocation(p1u,p1l,nump1,numt1,p1v)
            (p2loc,t2loc) = vectorLocation(p2u,p2l,nump2,numt2,p2v)
            loc12 = CartesianIndex(p1loc,t1loc,p2loc,t2loc)

            SAtally[loc12] += UInt32(1)
            SAtotal[loc12] += Sval
        
        end

    end

end # function