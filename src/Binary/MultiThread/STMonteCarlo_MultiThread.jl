
#= 
This module provides functions for MonteCarlo Integration of S and T Matrices
=#

"""
    STMonteCarlo_MultiThread!(Arrays,ArrayOfLocks,sigma,dsigmadt,UserInput)

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Brief
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3,4)
    - RandomPointMomentum for p ( species 1,2 only)
- Calculate T value 
- Take random points (u3,h3,p1,p2,u1,u2,h1,h2) and calculate valid p3 point/points 
- Calculate S3 value
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Take random points (u4,h3,p1,p2,u1,u2,h1,h2) and calculate valid p4 point/points 
- Calculate S4 value
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Update global S and T arrays with locks to prevent data races
"""
function STMonteCarlo_MultiThread!(GainTotal3::Array{Float64,9},GainTotal4::Array{Float64,9},LossTotal::Array{Float64,6},GainTally3::Array{UInt32,9},GainTally4::Array{UInt32,9},LossTally::Array{UInt32,6},ArrayOfLocks,sigma::Function,dsigmadt::Function,Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},numTiterPerThread::Int64,numSiterPerThread::Int64,scale::Float64,prog::Progress,thread_id::Int64)

    # Set Parameters

    (name1,name2,name3,name4,m1,m2,m3,m4,p1_low,p1_up,p1_grid_st,p1_num,u1_grid_st,u1_num,h1_grid_st,h1_num,p2_low,p2_up,p2_grid_st,p2_num,u2_grid_st,u2_num,h2_grid_st,h2_num,p3_low,p3_up,p3_grid_st,p3_num,u3_grid_st,u3_num,h3_grid_st,h3_num,p4_low,p4_up,p4_grid_st,p4_num,u4_grid_st,u4_num,h4_grid_st,h4_num) = Parameters
    #=Parameters = userInput.Parameters
    numTiterPerThread = userInput.numTiter
    numSiterPerThread = userInput.numSiter
    MinMax = userInput.MinMax

    name1 = Parameters.name1
    name2 = Parameters.name2
    name3 = Parameters.name3
    name4 = Parameters.name4
    m1 = Parameters.m1
    m2 = Parameters.m2
    m3 = Parameters.m3
    m4 = Parameters.m4

    p1_low = Parameters.p1_low
    p1_up = Parameters.p1_up
    p1_grid_st = Parameters.p1_grid
    p1_num = Parameters.p1_num

    u1_grid_st = Parameters.u1_grid
    u1_num = Parameters.u1_num

    h1_grid_st = Parameters.h1_grid
    h1_num = Parameters.h1_num

    p2_low = Parameters.p2_low
    p2_up = Parameters.p2_up
    p2_grid_st = Parameters.p2_grid
    p2_num = Parameters.p2_num

    u2_grid_st = Parameters.u2_grid
    u2_num = Parameters.u2_num

    h2_grid_st = Parameters.h2_grid
    h2_num = Parameters.h2_num

    p3_low = Parameters.p3_low
    p3_up = Parameters.p3_up
    p3_grid_st = Parameters.p3_grid
    p3_num = Parameters.p3_num

    u3_grid_st = Parameters.u3_grid
    u3_num = Parameters.u3_num

    h3_grid_st = Parameters.h3_grid
    h3_num = Parameters.h3_num

    p4_low = Parameters.p4_low
    p4_up = Parameters.p4_up
    p4_grid_st = Parameters.p4_grid
    p4_num = Parameters.p4_num

    u4_grid_st = Parameters.u4_grid
    u4_num = Parameters.u4_num

    h4_grid_st = Parameters.h4_grid
    h4_num = Parameters.h4_num

    # Set Arrays
    GainTotal3 = Arrays.GainTotal3
    GainTotal4 = Arrays.GainTotal4
    TAtotal = Arrays.TAtotal
    GainTally3 = Arrays.GainTally3
    GainTally4 = Arrays.GainTally4
    TAtally = Arrays.TAtally
    if MinMax
        p3Max = Arrays.p3Max
        u3MinMax = Arrays.u3MinMax
        p4Max = Arrays.p4Max
        u4MinMax = Arrays.u4MinMax
    end
    =#

    # Set up worker
    Threads.@spawn begin

    # allocate arrays for each thread
    p1v::Vector{Float64} = zeros(Float64,4)
    p1cv::Vector{Float64} = zeros(Float64,4)
    p1cBv::Vector{Float64} = zeros(Float64,4)
    p2v::Vector{Float64} = zeros(Float64,4)
    #pv::Vector{Float64} = zeros(Float64,3)
    pv::Vector{Float64} = zeros(Float64,4)
    βv::Vector{Float64} = zeros(Float64,3)
    #p3v::Vector{Float64} = zeros(Float64,3)
    p3v::Vector{Float64} = zeros(Float64,4)
    #p3pv::Vector{Float64} = zeros(Float64,3)
    p3pv::Vector{Float64} = zeros(Float64,4)
    #p4v::Vector{Float64} = zeros(Float64,3)
    p4v::Vector{Float64} = zeros(Float64,4)
    #p4pv::Vector{Float64} = zeros(Float64,3)
    p4pv::Vector{Float64} = zeros(Float64,4)
    Sval::Float64 = 0e0
    Svalp::Float64 = 0e0
    Tval::Float64 = 0e0
    p3_physical::Bool = true
    p3p_physical::Bool = true
    p4_physical::Bool = true
    p4p_physical::Bool = true
    NumStates::Int64 = 2

    p1loc::Int64 = 0
    p2loc::Int64 = 0
    u1loc::Int64 = 0
    u2loc::Int64 = 0
    h1loc::Int64 = 0
    h2loc::Int64 = 0
    p3loc::Int64 = 0
    u3loc::Int64 = 0
    u3locMirror::Int64 = 0
    h3loc::Int64 = 0
    h3locMirror::Int64 = 0
    p3ploc::Int64 = 0
    u3ploc::Int64 = 0
    h3ploc::Int64 = 0
    p4loc::Int64 = 0
    u4loc::Int64 = 0
    u4locMirror::Int64 = 0
    h4loc::Int64 = 0
    h4locMirror::Int64 = 0
    p4ploc::Int64 = 0
    u4ploc::Int64 = 0
    h4ploc::Int64 = 0
    loc12::CartesianIndex{6} = CartesianIndex(0,0,0,0,0,0)
    WeightFactors::Tuple{Float64,Float64,Float64,Float64} = (0e0,0e0,0e0,0e0)

    p1_grid::GridType = Grid_String_to_Type(p1_grid_st)
    p2_grid::GridType = Grid_String_to_Type(p2_grid_st)
    p3_grid::GridType = Grid_String_to_Type(p3_grid_st)
    p4_grid::GridType = Grid_String_to_Type(p4_grid_st)
    u1_grid::GridType = Grid_String_to_Type(u1_grid_st)
    u2_grid::GridType = Grid_String_to_Type(u2_grid_st)
    u3_grid::GridType = Grid_String_to_Type(u3_grid_st)
    u4_grid::GridType = Grid_String_to_Type(u4_grid_st)
    h1_grid::GridType = Grid_String_to_Type(h1_grid_st)
    h2_grid::GridType = Grid_String_to_Type(h2_grid_st)
    h3_grid::GridType = Grid_String_to_Type(h3_grid_st)
    h4_grid::GridType = Grid_String_to_Type(h4_grid_st)

    SmallParameters = (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,m1,m2,m3,m4)

    localGainTotal3::Array{Float64,3} = zeros(Float64,size(GainTotal3)[1:3])
    localGainTally3::Array{UInt32,3} = zeros(UInt32,size(GainTally3)[1:3])
    localGainTotal4::Array{Float64,3} = zeros(Float64,size(GainTotal4)[1:3])
    localGainTally4::Array{UInt32,3} = zeros(UInt32,size(GainTally4)[1:3])

    @inbounds for _ in 1:numTiterPerThread
        
        # generate p1 and p2 vectors initially as to not have to re-calculate, but not p2 magnitude as we need one free parameter to vary
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)

        RPointLogMomentum!(p1v,p1_up,p1_low,p1_num)
        RPointLogMomentum!(p2v,p2_up,p2_low,p2_num)

        # Tval
        (Tval,sBig,sSmol) = TValue(p1v,p2v,sigma,m1,m2,m3,m4)
        # Calculate T Array Location
        p1loc = location(p1_low,p1_up,p1_num,p1v[1],p1_grid)
        p2loc = location(p2_low,p2_up,p2_num,p2v[1],p2_grid)
        u1loc = location(u_low,u_up,u1_num,p1v[2],u1_grid)
        u2loc = location(u_low,u_up,u2_num,p2v[2],u2_grid)
        h1loc = location(h_low,h_up,h1_num,p1v[3],h1_grid)
        h2loc = location(h_low,h_up,h2_num,p2v[3],h2_grid)
        loc12 = CartesianIndex(p1loc,u1loc,h1loc,p2loc,u2loc,h2loc)

        fill!(localGainTally3,UInt32(0))
        fill!(localGainTally4,UInt32(0))

        if Tval != 0e0 # i.e. it is a valid interaction state

            #WeightFactors = WeightedFactors(p1v,p2v,m1,m2,sBig,sSmol,scale) 
            WeightFactors = WeightedFactors2(p1v,p2v,m1,m2,m3,m4,sBig,sSmol,scale) 

            fill!(localGainTotal3,Float64(0))
            fill!(localGainTotal4,Float64(0))
                    
            @inbounds for _ in 1:numSiterPerThread


                ImportanceSampling4!(p1v,p2v,p3v,p4v,p3pv,p4pv,sBig,sSmol,dsigmadt,localGainTotal3,localGainTotal4,localGainTally3,localGainTally4,SmallParameters,WeightFactors)

            
            #=
                # generate random p direction for use in both p3 and p4 calculations
                RPointSphereCosThetaPhi!(pv)

            # ========= p3 ========= #

                #set random p3 direction 
                p3v .= pv
                p3pv .= pv

                # Calculate p3 value with checks
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,m1,m2,m3,m4)

                # S Array Tallies
                # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states.
                #if NumStates != 0
                u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                u3locMirror = location(u_low,u_up,u3_num,-p3v[2],u3_grid)
                h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
                h3locMirror = location(h_low,h_up,h3_num,mod(p3v[3]+1e0,2e0),h3_grid)
                localGainTally3[u3loc,h3loc] += UInt32(1)
                localGainTally3[u3locMirror,h3locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,m1,m2,m3,m4)
                        localGainTotal3[p3loc,u3loc,h3loc] += Sval
                        if MinMax
                            localp3Max[u3loc,h3loc] = max(localp3Max[u3loc,h3loc],p3v[1])
                            localu3Min[p3loc,h3loc] = min(localu3Min[p3loc,h3loc],p3v[2])
                            localu3Max[p3loc,h3loc] = max(localu3Max[p3loc,h3loc],p3v[2])
                        end
                    end
                end

                if NumStates == 2
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,m1,m2,m3,m4)
                        localGainTotal3[p3loc,u3loc,h3loc] += Sval
                        if MinMax
                            localp3Max[u3loc,h3loc] = max(localp3Max[u3loc,h3loc],p3v[1])
                            localu3Min[p3loc,h3loc] = min(localu3Min[p3loc,h3loc],p3v[2])
                            localu3Max[p3loc,h3loc] = max(localu3Max[p3loc,h3loc],p3v[2])
                        end
                    end
                    if p3p_physical
                        u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                        h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
                        p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                        Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,m1,m2,m3,m4)
                        localGainTotal3[p3ploc,u3ploc,h3ploc] += Svalp
                        if MinMax
                            localp3Max[u3ploc,h3ploc] = max(localp3Max[u3ploc,h3ploc],p3pv[1])
                            localu3Min[p3ploc,h3ploc] = min(localu3Min[p3ploc,h3ploc],p3pv[2])
                            localu3Max[p3ploc,h3ploc] = max(localu3Max[p3ploc,h3ploc],p3pv[2])
                        end
                    end
                end

            # ========= p4 ========= #
                # only need to swap m3 for m4

                #set random p4 direction 
                p4v .= pv
                p4pv .= pv

                # Calculate p3 value with checks
                (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p2v,p1v,m2,m1,m4,m3)

                # S Array Tallies
                # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states.
                u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
                u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
                h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
                localGainTally4[u4loc,h4loc] += UInt32(1)
                localGainTally4[u4locMirror,h4locMirror] += UInt32(1)

                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,m1,m2,m3,m4)
                        localGainTotal4[p4loc,u4loc,h4loc] += Sval
                        if MinMax
                            localp4Max[u4loc,h4loc] = max(localp4Max[u4loc,h4loc],p4v[1])
                            localu4Min[p4loc,h4loc] = min(localu4Min[p4loc,h4loc],p4v[2])
                            localu4Max[p4loc,h4loc] = max(localu4Max[p4loc,h4loc],p4v[2])
                        end
                    end
                end

                if NumStates == 2
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,m1,m2,m3,m4)
                        localGainTotal4[p4loc,u4loc,h4loc] += Sval
                        if MinMax
                            localp4Max[u4loc,h4loc] = max(localp4Max[u4loc,h4loc],p4v[1])
                            localu4Min[p4loc,h4loc] = min(localu4Min[p4loc,h4loc],p4v[2])
                            localu4Max[p4loc,h4loc] = max(localu4Max[p4loc,h4loc],p4v[2])
                        end
                    end
                    if p4p_physical
                        u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                        h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
                        p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                        Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,m1,m2,m3,m4)
                        localGainTotal4[p4ploc,u4ploc,h4ploc] += Svalp
                        if MinMax
                            localp4Max[u4ploc,h4ploc] = max(localp4Max[u4ploc,h4ploc],p4pv[1])
                            localu4Min[p4ploc,h4ploc] = min(localu4Min[p4ploc,h4ploc],p4pv[2])
                            localu4Max[p4ploc,h4ploc] = max(localu4Max[p4ploc,h4ploc],p4pv[2])
                        end
                    end
                end
                =#

            end # S loop

        else # no valid interaction state
            # add one to tally of all relevant S tallies i.e. all momenta and all angles as no emission states are possible
            @view(localGainTally3[end,:,:]) .+= UInt32(1)
            @view(localGainTally4[end,:,:]) .+= UInt32(1)
        end

        # assign values to arrays
        @lock ArrayOfLocks[p1loc] begin
            LossTotal[loc12] += Tval
            LossTally[loc12] += UInt32(1)
            @view(GainTally3[:,:,:,loc12]) .+= localGainTally3
            @view(GainTally4[:,:,:,loc12]) .+= localGainTally4
            if Tval != 0e0
                @view(GainTotal3[:,:,:,loc12]) .+= localGainTotal3
                @view(GainTotal4[:,:,:,loc12]) .+= localGainTotal4
            end
        end

        if thread_id == 1 # on main thread
            next!(prog)
        end

    end # T loop

    end # Thread spawn 

end # function 
