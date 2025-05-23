"""
    SyncMonteCarloAxi_MultiThread!(SAtotal,SAtally,pMax,tMinMax,Parameters,numSiter)

# Arguments
- `SAtotal::Array{Float64,4}` : Array of stored integration totals for S matrix for 2+B->2+B+1 interaction
- `SAtally::Array{UInt32,4}` : Array of stored integration tallies for S matrix for 2+B->2+B+1 interaction
- `pMax::Array{Float64,3}` : Array of maximum momentum values for species 2
- `tMinMax::Array{Float64,3}` : Array of minimum and maximum theta values for species 2
- `Parameters::Tuple{Float64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numSiter::Int64` : Number of S iterations

# Output:
- Argument arrays SAtotal,SAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Brief
- Random Sample points in each of these domains
    - RandomPointSphere for theta (for species 1,2)
    - RandomPointMomentum for p ( species 1,2)
- Take random points (p1,p2,t1,t2) and calculate Synchrotron emissivity
- Find position in S arrays and allocated tallies and totals accordingly.
"""
function SyncMonteCarloAxi_MultiThread!(SAtotal::Array{Float64,4},SAtally::Array{UInt32,4},#=pMax::Array{Float64,5},tMinMax::Array{Float64,6}=#ArrayOfLocks,Parameters::Tuple{String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64, Float64},numTiter::Int64,numSiter::Int64,nThreads::Int64,prog::Progress)

    # Set Parameters
    (name1,name2,mu1,mu2,z1,z2,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,BMag) = Parameters

    # Set up workers
    Threads.@spawn begin

    # allocate arrays
    p1v::Vector{Float64} = zeros(Float64,2)
    p2v::Vector{Float64} = zeros(Float64,2)
    Sval::Float64 = 0e0

    p1loc::Int64 = 0
    p2loc::Int64 = 0
    u1loc::Int64 = 0
    u2loc::Int64 = 0

    localStotal::Array{Float64,2} = zeros(Float64,p1_num,u1_num)
    localStally::Array{UInt32,2} = zeros(UInt32,p1_num,u1_num)

    for _ in 1:numTiter

        # generate p2v (emitting particle)
        RPointSphereCosTheta!(p2v)
        RPointLogMomentum!(p2v,p2_low,p2_up,p2_num)
        p2loc = location(p2_low,p2_up,p2_num,p2v[1],p2_grid)
        u2loc = location(u_low,u_up,u2_num,p2v[2],u2_grid)

        fill!(localStotal,Float64(0))
        fill!(localStally,UInt32(0))

        for _ in 1:numSiter

            # generate p1v (photon)
            RPointSphereCosTheta!(p1v)
            RPointLogMomentum!(p1v,p1_low,p1_up,p1_num)

            # calculate S value
            Sval = SyncKernel(p1v,p2v,mu2,z2,BMag)
            # find S array location 
            p1loc = location(p1_low,p1_up,p1_num,p1v[1],p1_grid)
            u1loc = location(u_low,u_up,u1_num,p1v[2],u1_grid)

            localStally[p1loc,u1loc] += UInt32(1)
            localStotal[p1loc,u1loc] += Sval
        
        end

        # assign values to arrays
        @lock ArrayOfLocks[p2loc] begin
            @view(SAtally[:,:,p2loc,u2loc]) .+= localStally
            @view(SAtotal[:,:,p2loc,u2loc]) .+= localStotal
        end

        next!(prog)

    end

    end # workers 

end # function