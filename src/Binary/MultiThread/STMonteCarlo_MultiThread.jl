
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
function STMonteCarlo_MultiThread!(Arrays::ScatteringArrays,ArrayOfLocks,sigma::Function,dsigmadt::Function,userInput::BinaryUserInput,p)

    # Set Parameters
    Parameters = userInput.Parameters
    numTiterPerThread = userInput.numTiter
    numSiterPerThread = userInput.numSiter
    MinMax = userInput.MinMax

    name1 = Parameters.name1
    name2 = Parameters.name2
    name3 = Parameters.name3
    name4 = Parameters.name4
    mu1 = Parameters.mu1
    mu2 = Parameters.mu2
    mu3 = Parameters.mu3
    mu4 = Parameters.mu4

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
    SAtotal3 = Arrays.SAtotal3
    SAtotal4 = Arrays.SAtotal4
    TAtotal = Arrays.TAtotal
    SAtally3 = Arrays.SAtally3
    SAtally4 = Arrays.SAtally4
    TAtally = Arrays.TAtally
    if MinMax
        p3Max = Arrays.p3Max
        u3MinMax = Arrays.u3MinMax
        p4Max = Arrays.p4Max
        u4MinMax = Arrays.u4MinMax
    end

    # Set up worker
    Threads.@spawn begin

    # allocate arrays for each thread
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

    p1_grid = Grid_String_to_Type(p1_grid_st)
    p2_grid = Grid_String_to_Type(p2_grid_st)
    p3_grid = Grid_String_to_Type(p3_grid_st)
    p4_grid = Grid_String_to_Type(p4_grid_st)
    u1_grid = Grid_String_to_Type(u1_grid_st)
    u2_grid = Grid_String_to_Type(u2_grid_st)
    u3_grid = Grid_String_to_Type(u3_grid_st)
    u4_grid = Grid_String_to_Type(u4_grid_st)
    h1_grid = Grid_String_to_Type(h1_grid_st)
    h2_grid = Grid_String_to_Type(h2_grid_st)
    h3_grid = Grid_String_to_Type(h3_grid_st)
    h4_grid = Grid_String_to_Type(h4_grid_st)

    localSAtotal3 = zeros(Float64,size(SAtotal3)[1:3])
    localSAtally3 = zeros(UInt32,size(SAtally3)[1:2])
    localSAtotal4 = zeros(Float64,size(SAtotal4)[1:3])
    localSAtally4 = zeros(UInt32,size(SAtally4)[1:2])
    if MinMax
        localp3Max = zeros(Float64,size(p3Max)[1:2])
        localu3Min = zeros(Float64,size(u3MinMax)[2:3])
        localu3Max = zeros(Float64,size(u3MinMax)[2:3])
        localp4Max = zeros(Float64,size(p4Max)[1:2])
        localu4Min = zeros(Float64,size(u4MinMax)[2:3])
        localu4Max = zeros(Float64,size(u4MinMax)[2:3])
    end

    for _ in 1:numTiterPerThread

        # generate p1 and p2 vectors initially as to not have to re-calculate, but not p2 magnitude as we need one free parameter to vary
        RPointSphereCosThetaPhi!(p1v)
        RPointSphereCosThetaPhi!(p2v)

        RPointLogMomentum!(p1v,p1_up,p1_low,p1_num)
        RPointLogMomentum!(p2v,p2_up,p2_low,p2_num)

        # Tval
        Tval = TValue(p1v,p2v,sigma,mu1,mu2,mu3,mu4)
        # Calculate T Array Location
        p1loc = location(p1_low,p1_up,p1_num,p1v[1],p1_grid)
        p2loc = location(p2_low,p2_up,p2_num,p2v[1],p2_grid)
        u1loc = location(u_low,u_up,u1_num,p1v[2],u1_grid)
        u2loc = location(u_low,u_up,u2_num,p2v[2],u2_grid)
        h1loc = location(h_low,h_up,h1_num,p1v[3],h1_grid)
        h2loc = location(h_low,h_up,h2_num,p2v[3],h2_grid)
        loc12 = CartesianIndex(p1loc,u1loc,h1loc,p2loc,u2loc,h2loc)

        fill!(localSAtally3,UInt32(0))
        fill!(localSAtally4,UInt32(0))

        if Tval != 0e0 # i.e. it is a valid interaction state

            fill!(localSAtotal3,Float64(0))
            fill!(localSAtotal4,Float64(0))
            if MinMax
                fill!(localp3Max,Float64(0))
                fill!(localu3Min,Float64(0))
                fill!(localu3Max,Float64(0))
                fill!(localp4Max,Float64(0))
                fill!(localu4Min,Float64(0))
                fill!(localu4Max,Float64(0))
            end
                    
            @inbounds for _ in 1:numSiterPerThread

                # generate random p direction for use in both p3 and p4 calculations
                RPointSphereCosThetaPhi!(pv)

            # ========= p3 ========= #

                #set random p3 direction 
                p3v .= pv
                p3pv .= pv

                # Calculate p3 value with checks
                (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

                # S Array Tallies
                # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states.
                #if NumStates != 0
                u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                u3locMirror = location(u_low,u_up,u3_num,-p3v[2],u3_grid)
                h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
                h3locMirror = location(h_low,h_up,h3_num,mod(p3v[3]+1e0,2e0),h3_grid)
                    localSAtally3[u3loc,h3loc] += UInt32(1)
                    localSAtally3[u3locMirror,h3locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3loc,u3loc,h3loc] += Sval
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
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3loc,u3loc,h3loc] += Sval
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
                        Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        if MinMax
                            localSAtotal3[p3ploc,u3ploc,h3ploc] += Svalp
                            localp3Max[u3ploc,h3ploc] = max(localp3Max[u3ploc,h3ploc],p3pv[1])
                            localu3Min[p3ploc,h3ploc] = min(localu3Min[p3ploc,h3ploc],p3pv[2])
                            localu3Max[p3ploc,h3ploc] = max(localu3Max[p3ploc,h3ploc],p3pv[2])
                        end
                    end
                end

            # ========= p4 ========= #
                # only need to swap mu3 for mu4

                #set random p4 direction 
                p4v .= pv
                p4pv .= pv

                # Calculate p3 value with checks
                (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p1v,p2v,mu1,mu2,mu4,mu3)

                # S Array Tallies
                # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states.
                u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
                u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
                h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
                localSAtally4[u4loc,h4loc] += UInt32(1)
                localSAtally4[u4locMirror,h4locMirror] += UInt32(1)

                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4loc,u4loc,h4loc] += Sval
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
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4loc,u4loc,h4loc] += Sval
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
                        Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4ploc,u4ploc,h4ploc] += Svalp
                        if MinMax
                            localp4Max[u4ploc,h4ploc] = max(localp4Max[u4ploc,h4ploc],p4pv[1])
                            localu4Min[p4ploc,h4ploc] = min(localu4Min[p4ploc,h4ploc],p4pv[2])
                            localu4Max[p4ploc,h4ploc] = max(localu4Max[p4ploc,h4ploc],p4pv[2])
                        end
                    end
                end

            end # S loop

        else # no valid interaction state
            # add one to tally of all relevant S tallies i.e. all momenta and all angles as no emission states are possible
            localSAtally3 .+= UInt32(1)
            localSAtally4 .+= UInt32(1)
        end

        # assign values to arrays
        @lock ArrayOfLocks[p1loc] begin
            TAtotal[loc12] += Tval
            TAtally[loc12] += UInt32(1)
            @view(SAtally3[:,:,loc12]) .+= localSAtally3
            @view(SAtally4[:,:,loc12]) .+= localSAtally4
            if Tval != 0e0
                @view(SAtotal3[:,:,:,loc12]) .+= localSAtotal3
                @view(SAtotal4[:,:,:,loc12]) .+= localSAtotal4
                if MinMax
                    @view(p3Max[:,:,loc12]) .= max.(@view(p3Max[:,:,loc12]),localp3Max)
                    @view(u3MinMax[1,:,:,loc12]) .= min.(@view(u3MinMax[1,:,:,loc12]),localu3Min)
                    @view(u3MinMax[2,:,:,loc12]) .= max.(@view(u3MinMax[2,:,:,loc12]),localu3Max)
                    @view(p4Max[:,:,loc12]) .= max.(@view(p4Max[:,:,loc12]),localp4Max)
                    @view(u4MinMax[1,:,:,loc12]) .= min.(@view(u4MinMax[1,:,:,loc12]),localu4Min)
                    @view(u4MinMax[2,:,:,loc12]) .= max.(@view(u4MinMax[2,:,:,loc12]),localu4Max)
                end
            end 
        end 

    end # T loop

    next!(p)

    end # Thread spawn 

end # function 
