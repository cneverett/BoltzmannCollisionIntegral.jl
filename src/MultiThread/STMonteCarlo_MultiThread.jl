
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

"""
    STMonteCarloAxi_MultiThread!(SAtotal3,SAtotal4,TAtotal,SAtally3,SAtally4,TAtally,p3Max,t3MinMax,p4Max,t4MinMax,sigma,dsigmadt,Parameters,numTiterPerThread,numSiterPerThread)

# Arguments
- `SAtotal3::Array{Float64,6}` : Array of stored integration totals for S matrix for 12->34 interaction
- `SAtotal4::Array{Float64,6}` : Array of stored integration totals for S matrix for 12->43 interaction
- `TAtotal::Array{Float64,4}` : Array of stored integration totals for T matrix
- `SAtally3::Array{UInt32,5}` : Array of stored integration tallies for S matrix for 12->34 interaction
- `SAtally4::Array{UInt32,5}` : Array of stored integration tallies for S matrix for 12->43 interaction
- `TAtally::Array{UInt32,4}` : Array of stored integration tallies for T matrix
- `p3Max::Array{Float64,5}` : Array of maximum momentum values for species 3
- `t3MinMax::Array{Float64,6}` : Array of minimum and maximum theta values for species 3
- `p4Max::Array{Float64,5}` : Array of maximum momentum values for species 4
- `t4MinMax::Array{Float64,6}` : Array of minimum and maximum theta values for species 4
- `sigma::Function` : Cross section function for the interaction
- `dsigmadt::Function` : Differential cross section function for the interaction
- `Parameters::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numTiterPerThread::Int64` : Number of T iterations per thread
- `numSiterPerThread::Int64` : Number of S iterations per thread

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Breif
- Set up worker threads
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3,4)
    - RandomPointMomentum for p ( species 1,2 only)
- Take random points (t3,h3,p1,p2,t1,t2,h1,h2) and calculate valid p3 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Take random points (t4,h3,p1,p2,t1,t2,h1,h2) and calculate valid p4 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Update global S and T arrays with locks to prevent data races
"""
function STMonteCarloAxi_MultiThread!(SAtotal3::Array{Float64,6},SAtotal4::Array{Float64,6},TAtotal::Array{Float64,4},SAtally3::Array{UInt32,5},SAtally4::Array{UInt32,5},TAtally::Array{UInt32,4},ArrayOfLocks,p3Max::Array{Float64,5},p4Max::Array{Float64,5},t3MinMax::Array{Float64,6},t4MinMax::Array{Float64,6},sigma::Function,dsigmadt::Function,Parameters::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64,Int64,Int64},numTiterPerThread::Int64,numSiterPerThread::Int64)

    # Set Parameters
    (mu1,mu2,mu3,mu4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2) = Parameters

    # Set up worker
    Threads.@spawn begin

    # allocate arrays for each thread
    p1v::Vector{Float64} = zeros(Float64,3)
    p2v::Vector{Float64} = zeros(Float64,3)
    p3v::Vector{Float64} = zeros(Float64,3)
    p3pv::Vector{Float64} = zeros(Float64,3)
    p4v::Vector{Float64} = zeros(Float64,3)
    p4pv::Vector{Float64} = zeros(Float64,3)
    Sval::Float64 = 0e0
    Svalp::Float64 = 0e0
    Tval::Float64 = 0e0
    p3_physical::Bool = true
    p3p_physical::Bool = true
    p4_physical::Bool = true
    p4p_physical::Bool = true
    NumStates::Int64 = 2

    localSAtotal3 = zeros(Float64,size(SAtotal3)[1:2])
    localSAtally3 = zeros(UInt32,size(SAtally3)[1])
    localSAtotal4 = zeros(Float64,size(SAtotal4)[1:2])
    localSAtally4 = zeros(UInt32,size(SAtally4)[1])
    localp3Max = zeros(Float64,size(p3Max)[1])
    localt3Min = zeros(Float64,size(t3MinMax)[2])
    localt3Max = zeros(Float64,size(t3MinMax)[2])
    localp4Max = zeros(Float64,size(p4Max)[1])
    localt4Min = zeros(Float64,size(t4MinMax)[2])
    localt4Max = zeros(Float64,size(t4MinMax)[2])

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

        fill!(localSAtally3,UInt32(0))
        fill!(localSAtally4,UInt32(0))

        if Tval != 0e0 # i.e. it is a valid interaction state

            fill!(localSAtotal3,0e0)
            fill!(localp3Max,Float64(0))
            fill!(localt3Min,Float64(0))
            fill!(localt3Max,Float64(0))
            fill!(localSAtotal4,0e0)
            fill!(localp4Max,Float64(0))
            fill!(localt4Min,Float64(0))
            fill!(localt4Max,Float64(0))
                    
            @inbounds for _ in 1:numSiterPerThread

            # ========= p3 ========= #

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
                    localSAtally3[t3loc] += UInt32(1)
                    localSAtally3[t3locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3loc,t3loc] += Sval
                        localp3Max[t3loc] = max(localp3Max[t3loc],p3v[1])
                        localt3Min[p3loc] = min(localt3Min[p3loc],p3v[2])
                        localt3Max[p3loc] = max(localt3Max[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    t3ploc = location_t(numt3,p3pv[2])
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3loc,t3loc] += Sval
                        localp3Max[t3loc] = max(localp3Max[t3loc],p3v[1])
                        localt3Min[p3loc] = min(localt3Min[p3loc],p3v[2])
                        localt3Max[p3loc] = max(localt3Max[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = location_p3(p3u,p3l,nump3,p3pv[1])
                        Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3ploc,t3ploc] += Svalp
                        localp3Max[t3ploc] = max(localp3Max[t3ploc],p3pv[1])
                        localt3Min[p3ploc] = min(localt3Min[p3ploc],p3pv[2])
                        localt3Max[p3ploc] = max(localt3Max[p3ploc],p3pv[2])
                    end
                end

            # ========= p4 ========= #
                # only need to swap mu3 for mu4

                #generate random p4 direction 
                RPointSphereCosThetaPhi!(p4v)
                p4pv .= p4v

                # Calculate p3 value with checks
                (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p1v,p2v,mu1,mu2,mu4,mu3)

                # S Array Tallies
                # For each t3 sampled, p3 will be + or -ve, corresponding to a change in sign of t3. Therefore by sampling one t3 we are actually sampling t3 and -t3 with one or both having valid p3 states.
                #if NumStates != 0
                    t4loc = location_t(numt4,p4v[2])
                    t4locMirror = location_t(numt4,-p4v[2])
                    localSAtally4[t4loc] += UInt32(1)
                    localSAtally4[t4locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location_p3(p4u,p4l,nump4,p4v[1])
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4loc,t4loc] += Sval
                        localp4Max[t4loc] = max(localp4Max[t4loc],p4v[1])
                        localt4Min[p4loc] = min(localt4Min[p4loc],p4v[2])
                        localt4Max[p4loc] = max(localt4Max[p4loc],p4v[2])
                    end
                end

                if NumStates == 2
                    t4ploc = location_t(numt4,p4pv[2])
                    if p4_physical
                        p4loc = location_p3(p4u,p4l,nump4,p4v[1])
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4loc,t4loc] += Sval
                        localp4Max[t4loc] = max(localp4Max[t4loc],p4v[1])
                        localt4Min[p4loc] = min(localt4Min[p4loc],p4v[2])
                        localt4Max[p4loc] = max(localt4Max[p4loc],p4v[2])
                    end
                    if p4p_physical
                        p4ploc = location_p3(p4u,p4l,nump4,p4pv[1])
                        Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4ploc,t4ploc] += Svalp
                        localp4Max[t4ploc] = max(localp4Max[t4ploc],p4pv[1])
                        localt4Min[p4ploc] = min(localt4Min[p4ploc],p4pv[2])
                        localt4Max[p4ploc] = max(localt4Max[p4ploc],p4pv[2])
                    end
                end

            end # Sloop

        else # no valid interaction state
            # add one to tally of all relavant S tallies i.e. all momenta and all angles as no emission states are possible
            localSAtally3 .+= UInt32(1)
            localSAtally4 .+= UInt32(1)
        end

        # assign values to arrays
        @lock ArrayOfLocks[p1loc] begin
            TAtotal[loc12] += Tval
            TAtally[loc12] += UInt32(1)
            @view(SAtotal3[:,:,loc12]) .+= localSAtotal3
            @view(SAtally3[:,loc12]) .+= localSAtally3
            @view(SAtotal4[:,:,loc12]) .+= localSAtotal4
            @view(SAtally4[:,loc12]) .+= localSAtally4
            if Tval != 0e0
                @view(p3Max[:,loc12]) .= max.(@view(p3Max[:,loc12]),localp3Max)
                @view(t3MinMax[1,:,loc12]) .= min.(@view(t3MinMax[1,:,loc12]),localt3Min)
                @view(t3MinMax[2,:,loc12]) .= max.(@view(t3MinMax[2,:,loc12]),localt3Max)
                @view(p4Max[:,loc12]) .= max.(@view(p4Max[:,loc12]),localp4Max)
                @view(t4MinMax[1,:,loc12]) .= min.(@view(t4MinMax[1,:,loc12]),localt4Min)
                @view(t4MinMax[2,:,loc12]) .= max.(@view(t4MinMax[2,:,loc12]),localt4Max)
            end 
        end 

    end # Tloop

    end # Thread spwan 

end # function 
