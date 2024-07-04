
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
        - S Array will have dimensions ((p+1)x2px3t) for axisymmetric
        - T Array will have dimensions (2px2t) for axisymmetric
            - the x2 counts for total and tally values
            - this is likely to take up large data storage
        - an alternative is impliment where arrays are reduced to 2D in size and S and tally arrays are generated as sparse array. The dimensions of the S Array are ((p+1)x(2p)x(3t))
            - there will be an extra 1 for overflow momenta i.e. array acts like [p3 i, p3 i+1, p3 i+2 .... p3 nump3, overflow]
            - points that land under the momentum domain i.e. underflow are assigned to the lowest momentum bin. This in effect "re-boosts" them up to the lowest momentum.
    
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

SAtotal = Array{Float32,6}(undef, nump3+1,numt3,nump1,numt1,nump2,numt2);
TAtotal = Array{Float32,4}(undef,nump1,numt1,nump2,numt2);
SAtally = Array{UInt32,5}(undef,numt3,nump1,numt1,nump2,numt2);
TAtally = Array{UInt32,4}(undef,nump1,numt1,nump2,numt2);
ArrayOfLocks = [Threads.SpinLock() for _ in 1:numt2] 

p1v = zeros(Float32,3)
p2v = zeros(Float32,3)
p3v = zeros(Float32,3,2)

using BenchmarkTools

@btime STMonteCarloAxi_MultiThread!(SAtotal,TAtotal,SAtally,TAtally,ArrayOfLocks)

# ------------- =#


#function STMonteCarloAxi_MultiThread!(SAtotal::Array{Float32,6},TAtotal::Array{Float32,4},SAtally::Array{UInt32,5},TAtally::Array{UInt32,4},ArrayOfLocks,p3Max::Array{Float32,5},t3MinMax::Array{Float32,6})

    # check arrays are correct size 
    #size(AStally) != ((nump3+1),numt3,nump1,numt1,nump2,numt2) && error("ASally Array improperly sized")
    #size(ATtally) != (nump1,numt1,nump2,numt2) && error("ATally Array improperly sized")
    #size(SAtotal) != ((nump3+1),numt3,nump1,numt1,nump2,numt2) && error("S Total Array improperly sized")
    #size(TAtotal) != (nump1,numt1,nump2,numt2) && error("tally Array improperly sized")

    # Set up worker

    #Threads.@spawn begin

    # allocate arrays for each thread
    p1v::Vector{Float32} = zeros(Float32,3)
    p2v::Vector{Float32} = zeros(Float32,3)
    p3::ComplexF32 = 0f0+0f0im
    p3p::ComplexF32 = 0f0+0f0im
    th3v::Vector{Float32} = zeros(Float32,2)
    th3pv::Vector{Float32} = zeros(Float32,2)
    Sval::Float32 = 0f0
    Svalp::Float32 = 0f0
    Tval::Float32 = 0f0
    sumTerms::Vector{Float32} = zeros(Float32,8)

    #localSAtotal = zeros(Float32,size(SAtotal)[1:2])
    #localSAtally = zeros(UInt32,size(SAtally)[1])
    #localp3Max = zeros(Float32,size(p3Max)[1])
    #localt3Min = zeros(Float32,size(t3MinMax)[2])
    #localt3Max = zeros(Float32,size(t3MinMax)[2])

    @btime for _ in 1:numTiterPerThread

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

        #fill!(localSAtally,UInt32(0))

        if Tval != 0f0 # i.e. it is a valid interaction state

            #fill!(localSAtotal,0f0)
            #fill!(localp3Max,Float32(0))
            #fill!(localt3Min,Float32(0))
            #fill!(localt3Max,Float32(0))
            
            @inbounds for _ in 1:numSiterPerThread

                #generate random p3 direction 
                RPointSphereCosThetaPhi_2Elem!(th3v)
                th3pv .= th3v

                # Calculate p3 value with checks
                (p3,p3p) = Momentum3Value3(th3v,p1v,p2v)

                # Calculate S Array Location
                if isreal(p3) # && isreal(p3p) only need to check if one is real
                    p3_re = real(p3)
                    p3p_re = real(p3p)
                    p12 = p1v[1]^2
                    p22 = p2v[1]^2
                    sqm1p1 = sqrt(p12+mu1^2)
                    sqm2p2 = sqrt(p22+mu2^2)

                    if p3_re < 0f0
                        # adjust th3 values
                        p3_re *= -1
                        th3v[1] *= -1
                        th3v[2] = mod(th3v[2]+1f0,2f0)
                        t3loc = location(t3u,t3l,numt3,th3v[1])
                        # add one to tally 
                        #localSAtally[t3loc] += UInt32(1)
                        # check if p3 is physical, if so add to total
                        if (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3_re^2/(sqrt(mu3^2+p3_re^2)+mu3) > mu3-mu2-mu1)
                            p3loc = locationp3(p3u,p3l,nump3,p3_re)
                            Sval = SValue2(p3_re,th3v,p1v,p2v,sumTerms)
                            #localSAtotal[p3loc,t3loc] += Sval
                        end
                        
                        if p3p_re < 0f0 
                            # th3p values same as th3, no need to adjust instead use th3
                            p3p_re *= -1
                            # check if p3p is physical and not identical to p3
                            if p3p_re != p3_re && (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3p_re^2/(sqrt(mu3^2+p3p_re^2)+mu3) > mu3-mu2-mu1)
                                p3ploc = locationp3(p3u,p3l,nump3,p3p_re)
                                Svalp = SValue2(p3p_re,th3v,p1v,p2v,sumTerms)
                                #localSAtotal[p3ploc,t3loc] += Svalp
                            end
                        end
                        if p3p_re > 0 
                            # no need to adjust th3p values, just need to calculate t3ploc
                            t3ploc = location(t3u,t3l,numt3,th3pv[1])
                            # add one to tally 
                            #localSAtally[t3ploc] += UInt32(1)
                            # check if p3p is physical
                            if (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3p_re^2/(sqrt(mu3^2+p3p_re^2)+mu3) > mu3-mu2-mu1)
                                p3ploc = locationp3(p3u,p3l,nump3,p3p_re)
                                Svalp = SValue2(p3p_re,th3pv,p1v,p2v,sumTerms)
                                #localSAtotal[p3ploc,t3ploc] += Svalp
                            end
                        end

                    elseif  p3_re > 0f0
                        # no need to adjuct th3 values
                        t3loc = location(t3u,t3l,numt3,th3v[1])
                        # add on to tally
                        #localSAtally[t3loc] += UInt32(1)
                        # if physical
                        if (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3_re^2/(sqrt(mu3^2+p3_re^2)+mu3) > mu3-mu2-mu1)
                            p3loc = locationp3(p3u,p3l,nump3,p3_re)
                            Sval = SValue2(p3_re,th3v,p1v,p2v,sumTerms)
                            #localSAtotal[p3loc,t3loc] += Sval
                        end

                        if p3p_re > 0f0 
                            # th3p values same as th3, do not add to tally
                            # check if p3p is physical and not identical to p3
                            if p3p_re != p3_re && (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3p_re^2/(sqrt(mu3^2+p3p_re^2)+mu3) > mu3-mu2-mu1)
                                p3ploc = locationp3(p3u,p3l,nump3,p3p_re)
                                Svalp = SValue2(p3p_re,th3v,p1v,p2v,sumTerms)
                                #localSAtotal[p3ploc,t3loc] += Svalp
                            end
                        end
                        if p3p_re < 0 
                            # adjust th3p values
                            p3p_re *= -1
                            th3pv[1] *= -1
                            th3pv[2] = mod(th3pv[2]+1f0,2f0)
                            t3ploc = location(t3u,t3l,numt3,th3pv[1])
                            # add one to tally 
                            #localSAtally[t3ploc] += UInt32(1)
                            # check if p3p is physical
                            if (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3p_re^2/(sqrt(mu3^2+p3p_re^2)+mu3) > mu3-mu2-mu1)
                                p3ploc = locationp3(p3u,p3l,nump3,p3p_re)
                                Svalp = SValue2(p3p_re,th3pv,p1v,p2v,sumTerms)
                                #localSAtotal[p3ploc,t3ploc] += Svalp
                            end
                        end

                    else # p3_re = 0f0
                        if p3p_re > 0f0 
                            t3ploc = location(t3u,t3l,numt3,th3pv[1])
                            # add 1 to tally
                            #localSAtally[t3ploc] += UInt32(1)
                            # check if p3p is physical and not identical to p3
                            if (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3p_re^2/(sqrt(mu3^2+p3p_re^2)+mu3) > mu3-mu2-mu1)
                                p3ploc = locationp3(p3u,p3l,nump3,p3p_re)
                                Svalp = SValue2(p3p_re,th3v,p1v,p2v,sumTerms)
                                #localSAtotal[p3ploc,t3loc] += Svalp
                            end
                        end
                        if p3p_re < 0 
                            # adjust th3p values
                            p3p_re *= -1
                            th3pv[1] *= -1
                            th3pv[2] = mod(th3pv[2]+1f0,2f0)
                            t3ploc = location(t3u,t3l,numt3,th3pv[1])
                            # add one to tally 
                            #localSAtally[t3ploc] += UInt32(1)
                            # check if p3p is physical
                            if (p12/(sqm1p1+mu1)+p22/(sqm2p2+mu2)-p3p_re^2/(sqrt(mu3^2+p3p_re^2)+mu3) > mu3-mu2-mu1)
                                p3ploc = locationp3(p3u,p3l,nump3,p3p_re)
                                Svalp = SValue2(p3p_re,th3pv,p1v,p2v,sumTerms)
                                #localSAtotal[p3ploc,t3ploc] += Svalp
                            end
                        end     
                    end

                elseif (isreal(p3)==false) # && (isreal(p3p)==false) # only need to check one
                    # if both complex then t3 must be the same for both
                    p3_re = real(p3)
                    if p3_re < 0f0
                        th3v[1] *= -1
                    end
                    t3loc = location(t3u,t3l,numt3,th3v[1])
                    #localSAtally[t3loc] += UInt32(1)
                end

            end # Sloop

        else # no valid interaction state
            # add one to tally of all relavant S tallies i.e. all momenta and all angles as no emission states are possible
            #localSAtally .+= UInt32(1)
        end

        # assign values to arrays
        #=@lock ArrayOfLocks[p1loc] begin
            TAtotal[loc12] += Tval
            TAtally[loc12] += UInt32(1)
            @view(SAtotal[:,:,loc12]) .+= localSAtotal
            @view(SAtally[:,loc12]) .+= localSAtally
            if Tval != 0f0
                @view(p3Max[:,loc12]) .= max.(@view(p3Max[:,loc12]),localp3Max)
                @view(t3MinMax[1,:,loc12]) .= min.(@view(t3MinMax[1,:,loc12]),localt3Min)
                @view(t3MinMax[2,:,loc12]) .= max.(@view(t3MinMax[2,:,loc12]),localt3Max)
            end 
        end=# 

    end # Tloop

#    end # Thread spwan 

#end # function 