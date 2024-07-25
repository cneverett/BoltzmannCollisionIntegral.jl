
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

"""
    STMonteCarloAxi_Serial!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax,sigma,dsigmadt,Parameters,numTiter,numSiter)

# Arguments
- `SAtotal3::Array{Float32,6}` : Array of stored integration totals for S matrix for 12->34 interaction
- `SAtotal4::Array{Float32,6}` : Array of stored integration totals for S matrix for 12->43 interaction
- `TAtotal::Array{Float32,4}` : Array of stored integration totals for T matrix
- `SAtally3::Array{UInt32,5}` : Array of stored integration tallies for S matrix for 12->34 interaction
- `SAtally4::Array{UInt32,5}` : Array of stored integration tallies for S matrix for 12->43 interaction
- `TAtally::Array{UInt32,4}` : Array of stored integration tallies for T matrix
- `p3Max::Array{Float32,5}` : Array of maximum momentum values for species 3
- `t3MinMax::Array{Float32,6}` : Array of minimum and maximum theta values for species 3
- `p4Max::Array{Float32,5}` : Array of maximum momentum values for species 4
- `t4MinMax::Array{Float32,6}` : Array of minimum and maximum theta values for species 4
- `sigma::Function` : Cross section function for the interaction
- `dsigmadt::Function` : Differential cross section function for the interaction
- `Parameters::Tuple{Float32,Float32,Float32,Float32,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Int64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numTiter::Int64` : Number of T iterations
- `numSiter::Int64` : Number of S iterations

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Breif
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3,4)
    - RandomPointMomentum for p ( species 1,2 only)
- Take random points (t3,h3,p1,p2,t1,t2,h1,h2) and calculate valid p3 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Take random points (t4,h3,p1,p2,t1,t2,h1,h2) and calculate valid p4 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
"""
function STMonteCarloAxi_Serial!(SAtotal3::Array{Float32,6},SAtotal4::Array{Float32,6},TAtotal::Array{Float32,4},SAtally3::Array{UInt32,5},SAtally4::Array{UInt32,5},TAtally::Array{UInt32,4},p3Max::Array{Float32,5},p4Max::Array{Float32,5},t3MinMax::Array{Float32,6},t4MinMax::Array{Float32,6},sigma::Function,dsigmadt::Function,Parameters::Tuple{Float32,Float32,Float32,Float32,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Int64,Int64,Int64,Int64},numTiter::Int64,numSiter::Int64)

    # Set Parameters
    (mu1,mu2,mu3,mu4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2) = Parameters

    # allocate arrays
    p1v::Vector{Float32} = zeros(Float32,3)
    p2v::Vector{Float32} = zeros(Float32,3)
    p3v::Vector{Float32} = zeros(Float32,3)
    p3pv::Vector{Float32} = zeros(Float32,3)
    p4v::Vector{Float32} = zeros(Float32,3)
    p4pv::Vector{Float32} = zeros(Float32,3)
    Sval::Float32 = 0f0
    Svalp::Float32 = 0f0
    Tval::Float32 = 0f0
    p3_physical::Bool = true
    p3p_physical::Bool = true
    p4_physical::Bool = true
    p4p_physical::Bool = true
    NumStates::Int64 = 2

    for _ in 1:numTiter

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

        SAtotalView3 = @view SAtotal3[:,:,loc12]
        SAtotalView4 = @view SAtotal4[:,:,loc12]
        SAtallyView3 = @view SAtally3[:,loc12]
        SAtallyView4 = @view SAtally4[:,loc12]
        p3MaxView = @view p3Max[:,loc12]
        t3MinView = @view t3MinMax[1,:,loc12]
        t3MaxView = @view t3MinMax[2,:,loc12]
        p4MaxView = @view p4Max[:,loc12]
        t4MinView = @view t4MinMax[1,:,loc12]
        t4MaxView = @view t4MinMax[2,:,loc12]
        
        if Tval != 0f0 # i.e. it is a valid interaction state

            for _ in 1:numSiter # loop over a number of p3 orientations for a given p1 p2 state

            # === p3 === #
                #generate random p3 direction 
                RPointSphereCosThetaPhi!(p3v)
                p3pv .= p3v

                # Calculate p3 value
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

                # S Array Tallies
                # For each t3 sampled, p3 will be + or -ve, corresponding to a change in sign of t3. Therefore by sampling one t3 we are actually sampling t3 and -t3 with one or both having valid p3 states.
                #if NumStates != 0
                    t3loc = location_t(numt3,p3v[2])
                    t3locMirror = location_t(numt3,-p3v[2])
                    SAtallyView3[t3loc] += UInt32(1)
                    SAtallyView3[t3locMirror] += UInt32(1)
                #end
  
                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)
                        SAtotalView3[p3loc,t3loc] += Sval
                        p3MaxView[t3loc] = max(p3MaxView[t3loc],p3v[1])
                        t3MinView[p3loc] = min(t3MinView[p3loc],p3v[2])
                        t3MaxView[p3loc] = max(t3MaxView[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    t3ploc = location_t(numt3,p3pv[2])
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3)
                        SAtotalView3[p3loc,t3loc] += Sval
                        p3MaxView[t3loc] = max(p3MaxView[t3loc],p3v[1])
                        t3MinView[p3loc] = min(t3MinView[p3loc],p3v[2])
                        t3MaxView[p3loc] = max(t3MaxView[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = location_p3(p3u,p3l,nump3,p3pv[1])
                        Svalp = SValue(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3)
                        SAtotalView3[p3ploc,t3ploc] += Svalp
                        p3MaxView[t3ploc] = max(p3MaxView[t3ploc],p3pv[1])
                        t3MinView[p3ploc] = min(t3MinView[p3ploc],p3pv[2])
                        t3MaxView[p3ploc] = max(t3MaxView[p3ploc],p3pv[2])
                    end
                end

            # === p4 === #
                #generate random p4 direction 
                RPointSphereCosThetaPhi!(p4v)
                p4pv .= p4v

                # Calculate p4 value
                (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p1v,p2v,mu1,mu2,mu4,mu3)

                # S Array Tallies
                # For each t4 sampled, p4 will be + or -ve, corresponding to a change in sign of t4. Therefore by sampling one t4 we are actually sampling t4 and -t4 with one or both having valid p4 states.
                #if NumStates != 0
                    t4loc = location_t(numt4,p4v[2])
                    t4locMirror = location_t(numt4,-p4v[2])
                    SAtallyView4[t4loc] += UInt32(1)
                    SAtallyView4[t4locMirror] += UInt32(1)
                #end
    
                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location_p3(p4u,p4l,nump4,p4v[1])
                        Sval = SValue(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu4)
                        SAtotalView4[p4loc,t4loc] += Sval
                        p4MaxView[t4loc] = max(p4MaxView[t4loc],p4v[1])
                        t4MinView[p4loc] = min(t4MinView[p4loc],p4v[2])
                        t4MaxView[p4loc] = max(t4MaxView[p4loc],p4v[2])
                    end
                end

                if NumStates == 2
                    t4ploc = location_t(numt4,p4pv[2])
                    if p4_physical
                        p4loc = location_p3(p4u,p4l,nump4,p4v[1])
                        Sval = SValue(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu4)
                        SAtotalView4[p4loc,t4loc] += Sval
                        p4MaxView[t4loc] = max(p4MaxView[t4loc],p4v[1])
                        t4MinView[p4loc] = min(t4MinView[p4loc],p4v[2])
                        t4MaxView[p4loc] = max(t4MaxView[p4loc],p4v[2])
                    end
                    if p4p_physical
                        p4ploc = location_p3(p4u,p4l,nump4,p4pv[1])
                        Svalp = SValue(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu4)
                        SAtotalView4[p4ploc,t4ploc] += Svalp
                        p4MaxView[t4ploc] = max(p4MaxView[t4ploc],p4pv[1])
                        t4MinView[p4ploc] = min(t4MinView[p4ploc],p4pv[2])
                        t4MaxView[p4ploc] = max(t4MaxView[p4ploc],p4pv[2])
                    end
                end

            end # Sloop

        else # no valid interaction state
            # add one to tally of all relavant S tallies i.e. all momenta and all angles as no emission states are possible
            SAtallyView3 .+= UInt32(1)
            SAtallyView4 .+= UInt32(1)
        end

        # asign to T arrays
        TAtotal[loc12] += Tval # ST[3] doesn't change with S loop
        TAtally[loc12] += UInt32(1)

    end # Tloop

    return nothing

end