
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

"""
    STMonteCarloAxi_Serial!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax})

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

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Hidden Inputs (defined in Init.jl)
- Domain Boundaries (defined as CONST)
        - p bounds and divisions for species 1,3,4
        - theta divisions for species 1,3,4 ( bounds not needed as assumed [-1,1] )
        - phi divisions for species 1,3,4 ( bounds not needed as assumed [0,2] )
- Particle Masses (defined as CONST)
- numTiter and numSiter as the number of T and S integrations to perform.

# Calculation In Breif
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3)
    - RandomPointMomentum for p ( species 1,2 only )
- Take random points (t3,h1,p1,p2,t1,t2,h3,h4) and calculate valid p3 point/points 
- Find position in S and T arrays and allocated tallies and totals accordingly.
"""
function STMonteCarloAxi_Serial!(SAtotal::Array{Float32,6},TAtotal::Array{Float32,4},SAtally::Array{UInt32,5},TAtally::Array{UInt32,4},p3v::Vector{Float32},p3pv::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32},p3Max::Array{Float32,5},t3MinMax::Array{Float32,6},sigma::Function,dsigmadt::Function)


    for _ in 1:numTiter

        # generate p1 and p2 vectors initially as to not have to re-caculate, but not p2 magnitude as we need one free parameter to vary
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)

        RPointLogMomentum!(p1u,p1l,p1v,nump1)
        RPointLogMomentum!(p2u,p2l,p2v,nump2)
  
        # Tval
        Tval = TValue(p1v,p2v,sigma)
        # Calculate T Array Location
        (p1loc,t1loc) = vectorLocation(p1u,p1l,nump1,numt1,p1v)
        (p2loc,t2loc) = vectorLocation(p2u,p2l,nump2,numt2,p2v)
        loc12 = CartesianIndex(p1loc,t1loc,p2loc,t2loc)

        SAtotalView = @view SAtotal[:,:,loc12]
        SAtallyView = @view SAtally[:,loc12]
        p3MaxView = @view p3Max[:,loc12]
        t3MinView = @view t3MinMax[1,:,loc12]
        t3MaxView = @view t3MinMax[2,:,loc12]
        
        if Tval != 0f0 # i.e. it is a valid interaction state

            for _ in 1:numSiter # loop over a number of p3 orientations for a given p1 p2 state

                #generate random p3 direction 
                RPointSphereCosThetaPhi!(p3v)
                p3pv .= p3v

                # Calculate p3 value
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v)

                # S Array Tallies
                # For each t3 sampled, p3 will be + or -ve, corresponding to a change in sign of t3. Therefore by sampling one t3 we are actually sampling t3 and -t3 with one or both having valid p3 states.
                #if NumStates != 0
                    t3loc = location_t(numt3,p3v[2])
                    t3locMirror = location_t(numt3,-p3v[2])
                    SAtallyView[t3loc] += UInt32(1)
                    SAtallyView[t3locMirror] += UInt32(1)
                #end
  
                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v,dsigmadt)
                        SAtotalView[p3loc,t3loc] += Sval
                        p3MaxView[t3loc] = max(p3MaxView[t3loc],p3v[1])
                        t3MinView[p3loc] = min(t3MinView[p3loc],p3v[2])
                        t3MaxView[p3loc] = max(t3MaxView[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    t3ploc = location_t(numt3,p3pv[2])
                    if p3_physical
                        p3loc = location_p3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v,dsigmadt)
                        SAtotalView[p3loc,t3loc] += Sval
                        p3MaxView[t3loc] = max(p3MaxView[t3loc],p3v[1])
                        t3MinView[p3loc] = min(t3MinView[p3loc],p3v[2])
                        t3MaxView[p3loc] = max(t3MaxView[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = location_p3(p3u,p3l,nump3,p3pv[1])
                        Svalp = SValue(p3pv,p1v,p2v,dsigmadt)
                        SAtotalView[p3ploc,t3ploc] += Svalp
                        p3MaxView[t3ploc] = max(p3MaxView[t3ploc],p3pv[1])
                        t3MinView[p3ploc] = min(t3MinView[p3ploc],p3pv[2])
                        t3MaxView[p3ploc] = max(t3MaxView[p3ploc],p3pv[2])
                    end
                end

            end # Sloop

        else # no valid interaction state
            # add one to tally of all relavant S tallies i.e. all momenta and all angles as no emission states are possible
            SAtallyView .+= UInt32(1)
        end

        # asign to T arrays
        TAtotal[loc12] += Tval # ST[3] doesn't change with S loop
        TAtally[loc12] += UInt32(1)

    end # Tloop

    return nothing

end