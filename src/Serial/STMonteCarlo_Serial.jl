
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

#= 
Intput:
    - Domain Boundaries
        - p bounds and divisions for species 1,3,4
        - theta divisions for species 1,3,4 ( bounds not nessessarlity needed as assumed [0,1] )
        - phi divisions for species 1,3,4 ( bounds not nessessarlity needed as assumed [0,2] )

    - Particle Masses
        - for species 1,2,3,4

    - Array of stored integration totals and tallys
        - total is cumulative sum of reaction rate in that domain
        - tally is cumalitive total of points that have landed in that doimain
        - S Array will have dimensions (3px3t) for axisymmetric and (3px3tx3h) for anisotropic
        - T Array will have dimensions (2px2t) for axisymmetric
            - the x2 counts for total and tally values
            - this is likely to take up large data storage
        - an alternative is impliment where arrays are reduced to 2D in size and S and tally arrays are generated as sparse array. The number of elements is ((3p)x(3t))
            - T array not sparse as all entries should be non-zero
            - there will be an extra 2 entries one for underflow and one for overflow momenta i.e. array acts like [underflow, overflow, p3 i, p3 i+1, p3 i+2 .... p3 nump3]
    
    - n integration points

Calculation:
    - Random Sample points in each of these domains
        - RandomPointSphere for theta and phi
        - RandomPointMomentum for p ( species 3,4 only )
    - Take random points (t3,h1,p1,p2,t1,t2,h3,h4) and calculate valid p3 point/points 
        - if non-valid point at t3 h1 then this will add 0 to MC integration total and 1 to MC integration tally
        - if valid point calculate reaction rate value 
    - Find position in arry of stored values corresponding to integration point and add reaction rate value to integration total and 1 to integration tally
    - Run for n integration points

Output:
    - Array of stored integration totals and tallys
    - One array (S array) gives rate of reaction to particular state/particles 1(2) from state 34 i.e. rate of emission of 1 from reaction 34->1(2)
    = One array (T array) gives rate of reaction from state/particles 34 to any state 12 i.e. rate of absorption of 34 in reaction 34->12

=#

#= for testing --

# Dependancies
include("../Common/MyPhysicalConstants.jl")
include("../Common/ParticleData.jl")
include("../Common\\Init.jl")
include("../Common\\DifferentialCrossSectionFunctions.jl")
include("../Common\\Momentum3Values.jl")
include("../Common\\RandomPointMomentum.jl")
include("../Common\\RandomPointSphere.jl")
include("../Common/MandelstramChecks.jl")
include("../Common\\STValue.jl")
include("../Common/UsefulGridValueFunctions.jl")
include("../Common/PhaseSpaceFactors.jl")
include("../Common/Location.jl")

using BenchmarkTools
SAtotal = Array{Float32,6}(undef, nump3+1,numt3,nump1,numt1,nump2,numt2);
TAtotal = Array{Float32,4}(undef,nump1,numt1,nump2,numt2);
SAtally = Array{UInt32,5}(undef,numt3,nump1,numt1,nump2,numt2);
TAtally = Array{UInt32,4}(undef,nump1,numt1,nump2,numt2);
p1v= zeros(Float32,3)
p2v = zeros(Float32,3)
p3v = zeros(Float32,3,2)
p3vTest = zeros(Float32,3)
p3vpTest= zeros(Float32,3)
ST = zeros(Float32,3)
# -------------

@btime STMonteCarloAxi_Serial!(SAtotal,TAtotal,SAtally,TAtally,p3v,p1v,p2v,ST)

@btime STMonteCarloAxi_SerialTest!(SAtotal,TAtotal,SAtally,TAtally,p3v,p1v,p2v)
@btime STMonteCarloAxi_SerialTest!(SAtotal,TAtotal,SAtally,TAtally,p3vTest,p3vpTest,p1v,p2v)

=#


function STMonteCarloAxi_Serial!(SAtotal::Array{Float32,6},TAtotal::Array{Float32,4},SAtally::Array{UInt32,5},TAtally::Array{UInt32,4},p3v::Vector{Float32},p3pv::Vector{Float32},p1v::Vector{Float32},p2v::Vector{Float32},p3Max::Array{Float32,5},t3MinMax::Array{Float32,6})


    for _ in 1:numTiter

        # generate p1 and p2 vectors initially as to not have to re-caculate, but not p2 magnitude as we need one free parameter to vary
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)

        RPointLogMomentum!(p1u,p1l,p1v,nump1)
        RPointLogMomentum!(p2u,p2l,p2v,nump2)
  
        # Tval
        Tval = TValue(p1v,p2v)
        # Calculate T Array Location
        (p1loc,t1loc) = vectorLocation(p1u,p1l,t1u,t1l,nump1,numt1,p1v)
        (p2loc,t2loc) = vectorLocation(p2u,p2l,t2u,t2l,nump2,numt2,p2v)
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
                if NumStates != 0
                    t3loc = location(t3u,t3l,numt3,p3v[2])
                    t3locMirror = location(t3u,t3l,numt3,-p3v[2])
                    SAtallyView[t3loc] += UInt32(1)
                    SAtallyView[t3locMirror] += UInt32(1)
                end
  
                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = locationp3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v)
                        SAtotalView[p3loc,t3loc] += Sval
                        p3MaxView[t3loc] = max(p3MaxView[t3loc],p3v[1])
                        t3MinView[p3loc] = min(t3MinView[p3loc],p3v[2])
                        t3MaxView[p3loc] = max(t3MaxView[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    t3ploc = location(t3u,t3l,numt3,p3pv[2])
                    if p3_physical
                        p3loc = locationp3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v)
                        SAtotalView[p3loc,t3loc] += Sval
                        p3MaxView[t3loc] = max(p3MaxView[t3loc],p3v[1])
                        t3MinView[p3loc] = min(t3MinView[p3loc],p3v[2])
                        t3MaxView[p3loc] = max(t3MaxView[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = locationp3(p3u,p3l,nump3,p3pv[1])
                        Svalp = SValue(p3pv,p1v,p2v)
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