
#= 
This module provides functions for MonteCarlo Integration of S and T Matricies
=#

"""
    STMonteCarloAxi_MultiThread!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax})

Intput:

    - Domain Boundaries (defined as CONST in Init.jl)
        - p bounds and divisions for species 1,3,4
        - theta divisions for species 1,3,4 ( bounds not nessessarlity needed as assumed [0,1] )
        - phi divisions for species 1,3,4 ( bounds not nessessarlity needed as assumed [0,2] )
    - Particle Masses (defined as CONST in Init.jl)
        - for species 1,2,3,4
    - Array of stored integration totals and tallys 
        - total is cumulative sum of reaction rate in that domain
        - tally is cumalitive total of points that have been sampled in that doimain
        - S Array will have dimensions ((nump3+1) x numt3 x nump1 x numt1 x nump2 x numt2) for axisymmetric
            - extra entry for p3 is for overflow momenta i.e. array acts like [p3 i, p3 i+1, p3 i+2 .... p3 nump3, overflow]
        - T Array will have dimensions (nump1 x numt1 x nump2 x numt2) for axisymmetric
    - numTiter and numSiter (defined in Init.jl) as the number of T and S integrations to perform.
    - nThreads as the number of threads to run the integration on

Calculation:

    - Set up workers to perform the integration on multiple threads
    - Random Sample points in each of these domains
        - RandomPointSphere for theta and phi (for species 1,2,3)
        - RandomPointMomentum for p ( species 1,2 only )
    - Take random points (t3,h1,p1,p2,t1,t2,h3,h4) and calculate valid p3 point/points 
    - Find position in S and T arrays and allocated tallies and totals accordingly (using locks to ensure single thread access to arrays). 

Output:

    - Edited arrays of stored integration totals and tallys
    - One array (S array) gives rate of reaction to particular state/particles 1(2) from state 34 i.e. rate of emission of 1 from reaction 34->1(2)
    - One array (T array) gives rate of reaction from state/particles 34 to any state 12 i.e. rate of absorption of 34 in reaction 34->12

"""
function STMonteCarloAxi_MultiThread!(SAtotal::Array{Float32,6},TAtotal::Array{Float32,4},SAtally::Array{UInt32,5},TAtally::Array{UInt32,4},ArrayOfLocks,p3Max::Array{Float32,5},t3MinMax::Array{Float32,6})

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

        RPointLogMomentum!(p1u,p1l,p1v,nump1)
        RPointLogMomentum!(p2u,p2l,p2v,nump2)

        # Tval
        Tval = TValue(p1v,p2v)
        # Calculate T Array Location
        (p1loc,t1loc) = vectorLocation(p1u,p1l,t1u,t1l,nump1,numt1,p1v)
        (p2loc,t2loc) = vectorLocation(p2u,p2l,t2u,t2l,nump2,numt2,p2v)
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
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v)

                # S Array Tallies
                # For each t3 sampled, p3 will be + or -ve, corresponding to a change in sign of t3. Therefore by sampling one t3 we are actually sampling t3 and -t3 with one or both having valid p3 states.
                if NumStates != 0
                    t3loc = location(t3u,t3l,numt3,p3v[2])
                    t3locMirror = location(t3u,t3l,numt3,-p3v[2])
                    localSAtally[t3loc] += UInt32(1)
                    localSAtally[t3locMirror] += UInt32(1)
                end

                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = locationp3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v)
                        localSAtotal[p3loc,t3loc] += Sval
                        localp3Max[t3loc] = max(localp3Max[t3loc],p3v[1])
                        localt3Min[p3loc] = min(localt3Min[p3loc],p3v[2])
                        localt3Max[p3loc] = max(localt3Max[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    t3ploc = location(t3u,t3l,numt3,p3pv[2])
                    if p3_physical
                        p3loc = locationp3(p3u,p3l,nump3,p3v[1])
                        Sval = SValue(p3v,p1v,p2v)
                        localSAtotal[p3loc,t3loc] += Sval
                        localp3Max[t3loc] = max(localp3Max[t3loc],p3v[1])
                        localt3Min[p3loc] = min(localt3Min[p3loc],p3v[2])
                        localt3Max[p3loc] = max(localt3Max[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = locationp3(p3u,p3l,nump3,p3pv[1])
                        Svalp = SValue(p3pv,p1v,p2v)
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