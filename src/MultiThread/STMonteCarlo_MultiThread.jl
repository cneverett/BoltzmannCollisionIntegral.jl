
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

"""
    STMonteCarloAxi_MultiThread!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax,sigma,dsigmadt,Parameters,numTiterPerThread,numSiterPerThread)

# Arguments
- `SAtotal::Array{Float32,6}` : Array of stored integration totals for S matrix
- `TAtotal::Array{Float32,4}` : Array of stored integration totals for T matrix
- `SAtally::Array{UInt32,5}` : Array of stored integration tallies for S matrix
- `TAtally::Array{UInt32,4}` : Array of stored integration tallies for T matrix
- `p3v::Vector{Float32}` : Vector of momentum values for species 3
- `p3pv::Vector{Float32}` : Second Vector of momentum values for species 3 for when two states are possible
- `p1v::Vector{Float32}` : Vector of momentum values for species 1
- `p2v::Vector{Float32}` : Vector of momentum values for species 2
- `p3Max::Array{Float32,5}` : Array of maximum momentum values for species 3
- `t3MinMax::Array{Float32,6}` : Array of minimum and maximum theta values for species 3
- `sigma::Function` : Cross section function for the interaction
- `dsigmadt::Function` : Differential cross section function for the interaction
- `Parameters::Tuple{Float32,Float32,Float32,Float32,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numTiterPerThread::Int64` : Number of T iterations per thread
- `numSiterPerThread::Int64` : Number of S iterations per thread

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Hidden Inputs (defined in Init.jl)
- Domain Boundaries (defined as CONST)
        - p bounds and divisions for species 1,3,4
        - theta divisions for species 1,3,4 ( bounds not needed as assumed [-1,1] )
        - phi divisions for species 1,3,4 ( bounds not needed as assumed [0,2] )
- Particle Masses (defined as CONST)
- numTiterPerThread and numSiterPerThread as the number of T and S integrations to perform.

# Calculation In Breif
- Set up worker threads
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3)
    - RandomPointMomentum for p ( species 1,2 only )
- Take random points (t3,h1,p1,p2,t1,t2,h3,h4) and calculate valid p3 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Update global S and T arrays with locks to prevent data races
"""
function STMonteCarloAxi_MultiThread!(SAtotal::Array{Float32,6},TAtotal::Array{Float32,4},SAtally::Array{UInt32,5},TAtally::Array{UInt32,4},ArrayOfLocks,p3Max::Array{Float32,5},t3MinMax::Array{Float32,6},sigma::Function,dsigmadt::Function,Parameters::Tuple{Float32,Float32,Float32,Float32,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Int64,Int64,Int64},numTiterPerThread::Int64,numSiterPerThread::Int64)

    # Set Parameters
    (mu1,mu2,mu3,mu4,p3l,p3u,nump3,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt1,numt2) = Parameters

    # Set up worker
    Threads.@spawn begin

    # allocate arrays for each thread
    p1v::Vector{Float32} = zeros(Float32,3)
    p2v::Vector{Float32} = zeros(Float32,3)
    p3v::Vector{Float32} = zeros(Float32,3)
    p3pv::Vector{Float32} = zeros(Float32,3)
    Sval::Float32 = 0f0
    Svalp::Float32 = 0f0
    Tval::Float32 = 0f0
    p3_physical::Bool = true
    p3p_physical::Bool = true
    NumStates::Int64 = 2

    localSAtotal = zeros(Float32,size(SAtotal)[1:2])
    localSAtally = zeros(UInt32,size(SAtally)[1])
    localp3Max = zeros(Float32,size(p3Max)[1])
    localt3Min = zeros(Float32,size(t3MinMax)[2])
    localt3Max = zeros(Float32,size(t3MinMax)[2])

    for _ in 1:numTiterPerThread

        # generate p1 and p2 vectors initially as to not have to re-caculate, but not p2 magnitude as we need one free parameter to vary
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)

        RPointLogMomentum!(p1v,p1u,p1l,nump1)
        RPointLogMomentum!(p2v,p2u,p2l,nump2)

        # Tval
        Tval = TValue(p1v,p2v,sigma,mu1,mu2,mu3,mu4)
        # Calculate T Array Location
        (p1loc,t1loc) = vectorLocation(p1u,p1l,nump1,numt1,p1v)
        (p2loc,t2loc) = vectorLocation(p2u,p2l,nump2,numt2,p2v)
        loc12 = CartesianIndex(p1loc,t1loc,p2loc,t2loc)

        fill!(localSAtally,UInt32(0))

        if Tval != 0f0 # i.e. it is a valid interaction state

            fill!(localSAtotal,0f0)
            fill!(localp3Max,Float32(0))
            fill!(localt3Min,Float32(0))
            fill!(localt3Max,Float32(0))
                    
            @inbounds for _ in 1:numSiterPerThread

                #generate random p3 direction 
                RPointSphereCosThetaPhi!(p3v)
                p3pv .= p3v

                # Calculate p3 value with checks
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

                # S Array Tallies
                # For each t3 sampled, p3 will be + or -ve, corresponding to a change in sign of t3. Therefore by sampling one t3 we are actually sampling t3 and -t3 with one or both having valid p3 states.
                #if NumStates != 0
                    t3loc = location_t(numt3,p3v[2])
                    t3locMirror = location_t(numt3,-p3v[2])
                    localSAtally[t3loc] += UInt32(1)
                    localSAtally[t3locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)
                        localSAtotal[p3loc,t3loc] += Sval
                        localp3Max[t3loc] = max(localp3Max[t3loc],p3v[1])
                        localt3Min[p3loc] = min(localt3Min[p3loc],p3v[2])
                        localt3Max[p3loc] = max(localt3Max[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    t3ploc = location_t(numt3,p3pv[2])
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)
                        localSAtotal[p3loc,t3loc] += Sval
                        localp3Max[t3loc] = max(localp3Max[t3loc],p3v[1])
                        localt3Min[p3loc] = min(localt3Min[p3loc],p3v[2])
                        localt3Max[p3loc] = max(localt3Max[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = location_p3(p3u,p3l,nump3,p3pv[1])
                        Svalp = SValue(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3)
                        localSAtotal[p3ploc,t3ploc] += Svalp
                        localp3Max[t3ploc] = max(localp3Max[t3ploc],p3pv[1])
                        localt3Min[p3ploc] = min(localt3Min[p3ploc],p3pv[2])
                        localt3Max[p3ploc] = max(localt3Max[p3ploc],p3pv[2])
                    end
                end

            end # Sloop

        else # no valid interaction state
            # add one to tally of all relavant S tallies i.e. all momenta and all angles as no emission states are possible
            localSAtally .+= UInt32(1)
        end

        # assign values to arrays
        @lock ArrayOfLocks[p1loc] begin
            TAtotal[loc12] += Tval
            TAtally[loc12] += UInt32(1)
            @view(SAtotal[:,:,loc12]) .+= localSAtotal
            @view(SAtally[:,loc12]) .+= localSAtally
            if Tval != 0f0
                @view(p3Max[:,loc12]) .= max.(@view(p3Max[:,loc12]),localp3Max)
                @view(t3MinMax[1,:,loc12]) .= min.(@view(t3MinMax[1,:,loc12]),localt3Min)
                @view(t3MinMax[2,:,loc12]) .= max.(@view(t3MinMax[2,:,loc12]),localt3Max)
            end 
        end 

    end # Tloop

    end # Thread spwan 

end # function 
