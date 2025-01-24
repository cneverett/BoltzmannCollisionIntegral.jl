
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

"""
    STMonteCarloAxi_Serial!(SAtotal3,SAtotal4,TAtotal,SAtally3,SAtally4,TAtally,p3Max,p4Max,t3MinMax,t4MinMax,sigma,dsigmadt,Parameters,numTiter,numSiter)

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
- `numTiter::Int64` : Number of T iterations
- `numSiter::Int64` : Number of S iterations

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Breif
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3,4)
    - RandomPointMomentum for p ( species 1,2 only)
- Take random points (u3,h3,p1,p2,u1,u2,h1,h2) and calculate valid p3 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Take random points (u4,h3,p1,p2,u1,u2,h1,h2) and calculate valid p4 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
"""
function STMonteCarloAxi_Serial!(SAtotal3::Array{Float64,6},SAtotal4::Array{Float64,6},TAtotal::Array{Float64,4},SAtally3::Array{UInt32,5},SAtally4::Array{UInt32,5},TAtally::Array{UInt32,4},p3Max::Array{Float64,5},p4Max::Array{Float64,5},t3MinMax::Array{Float64,6},t4MinMax::Array{Float64,6},sigma::Function,dsigmadt::Function,Parameters::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64,Int64,Int64},numTiter::Int64,numSiter::Int64)

    # Set Parameters
    (mu1,mu2,mu3,mu4,p3_low,p3_up,num_p3,p4_low,p4_up,num_p4,p1_low,p1_up,num_p1,p2_low,p2_up,num_p2,num_u3,num_u4,num_u1,num_u2) = Parameters

    # allocate arrays
    p1v::Vector{Float64} = zeros(Float64,3)
    p2v::Vector{Float64} = zeros(Float64,3)
    pv::Vector{Float64} = zeros(Float64,3)
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

    for _ in 1:numTiter

        # generate p1 and p2 vectors initially as to not have to re-caculate
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)

        RPointLogMomentum!(p1v,p1_up,p1_low,num_p1)
        RPointLogMomentum!(p2v,p2_up,p2_low,num_p2)
  
        # Tval
        Tval = TValue(p1v,p2v,sigma,mu1,mu2,mu3,mu4)
        # Calculate T Array Location
        p1loc = location(p1_up,p1_low,num_p1,p1v[1],"l")
        p2loc = location(p2_up,p2_low,num_p2,p2v[1],"l")
        u1loc = location(u_up,u_low,num_u1,p1v[2],"u")
        u2loc = location(u_up,u_low,num_u2,p2v[2],"u")
        (p1loc,u1loc) = vectorLocation(p1_up,p1_low,num_p1,num_u1,p1v)
        (p2loc,u2loc) = vectorLocation(p2_up,p2_low,num_p2,num_u2,p2v)
        loc12 = CartesianIndex(p1loc,u1loc,p2loc,u2loc)

        SAtotalView3 = @view SAtotal3[:,:,loc12]
        SAtotalView4 = @view SAtotal4[:,:,loc12]
        SAtallyView3 = @view SAtally3[:,loc12]
        SAtallyView4 = @view SAtally4[:,loc12]
        p3MaxView = @view p3Max[:,loc12]
        u3MinView = @view t3MinMax[1,:,loc12]
        u3MaxView = @view t3MinMax[2,:,loc12]
        p4MaxView = @view p4Max[:,loc12]
        u4MinView = @view t4MinMax[1,:,loc12]
        u4MaxView = @view t4MinMax[2,:,loc12]
        
        if Tval != 0e0 # i.e. it is a valid interaction state

            for _ in 1:numSiter # loop over a number of p3 orientations for a given p1 p2 state

                # generate random p direction for use in both p3 and p4 calculations
                RPointSphereCosThetaPhi!(pv)

            # === p3 === #
                #set random p3 direction 
                p3v .= pv
                p3pv .= pv

                # Calculate p3 value
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

                # S Array Tallies
                # For each t3 sampled, p3 will be + or -ve, corresponding to a change in sign of t3. Therefore by sampling one t3 we are actually sampling t3 and -t3 with one or both having valid p3 states.
                #if NumStates != 0
                    u3loc = location(u_up,u_low,num_u3,p3v[2],"u")
                    u3locMirror = location(u_up,u_low,num_u3,-p3v[2],"u")
                    SAtallyView3[u3loc] += UInt32(1)
                    SAtallyView3[u3locMirror] += UInt32(1)
                #end
  
                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location(p3_up,p3_low,num_p3,p3v[1],"l")
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView3[p3loc,u3loc] += Sval
                        p3MaxView[u3loc] = max(p3MaxView[u3loc],p3v[1])
                        u3MinView[p3loc] = min(u3MinView[p3loc],p3v[2])
                        u3MaxView[p3loc] = max(u3MaxView[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    if p3_physical
                        p3loc = location(p3_up,p3_low,num_p3,p3v[1],"l")
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView3[p3loc,u3loc] += Sval
                        p3MaxView[u3loc] = max(p3MaxView[u3loc],p3v[1])
                        u3MinView[p3loc] = min(u3MinView[p3loc],p3v[2])
                        u3MaxView[p3loc] = max(u3MaxView[p3loc],p3v[2])
                    end
                    if p3p_physical
                        u3ploc = location(u_up,u_low,num_u3,p3pv[2],"u")
                        p3ploc = location(p3_up,p3_low,num_p3,p3pv[1],"l")
                        Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView3[p3ploc,u3ploc] += Svalp
                        p3MaxView[u3ploc] = max(p3MaxView[u3ploc],p3pv[1])
                        u3MinView[p3ploc] = min(u3MinView[p3ploc],p3pv[2])
                        u3MaxView[p3ploc] = max(u3MaxView[p3ploc],p3pv[2])
                    end
                end

            # === p4 === #
                #set random p4 direction 
                p4v .= pv
                p4pv .= pv

                # Calculate p4 value
                (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p2v,p1v,mu2,mu1,mu4,mu3)

                # S Array Tallies
                # For each t4 sampled, p4 will be + or -ve, corresponding to a change in sign of t4. Therefore by sampling one t4 we are actually sampling t4 and -t4 with one or both having valid p4 states.
                #if NumStates != 0
                    u4loc = location(u_up,u_low,num_u4,p4v[2],"u")
                    u4locMirror = location(u_up,u_low,num_u4,-p4v[2],"u")
                    SAtallyView4[u4loc] += UInt32(1)
                    SAtallyView4[u4locMirror] += UInt32(1)
                #end
    
                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location(p4_up,p4_low,num_p4,p4v[1],"l")
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu4,mu4)
                        SAtotalView4[p4loc,u4loc] += Sval
                        p4MaxView[u4loc] = max(p4MaxView[u4loc],p4v[1])
                        u4MinView[p4loc] = min(u4MinView[p4loc],p4v[2])
                        u4MaxView[p4loc] = max(u4MaxView[p4loc],p4v[2])
                    end
                end

                if NumStates == 2
                    if p4_physical
                        p4loc = location(p4_up,p4_low,num_p4,p4v[1],"l")
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView4[p4loc,u4loc] += Sval
                        p4MaxView[u4loc] = max(p4MaxView[u4loc],p4v[1])
                        u4MinView[p4loc] = min(u4MinView[p4loc],p4v[2])
                        u4MaxView[p4loc] = max(u4MaxView[p4loc],p4v[2])
                    end
                    if p4p_physical
                        t4ploc = location(u_up,u_low,num_u4,p4pv[2],"u")
                        p4ploc = location(p4_up,p4_low,num_p4,p4pv[1],"l")
                        Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView4[p4ploc,t4ploc] += Svalp
                        p4MaxView[t4ploc] = max(p4MaxView[t4ploc],p4pv[1])
                        u4MinView[p4ploc] = min(u4MinView[p4ploc],p4pv[2])
                        u4MaxView[p4ploc] = max(u4MaxView[p4ploc],p4pv[2])
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