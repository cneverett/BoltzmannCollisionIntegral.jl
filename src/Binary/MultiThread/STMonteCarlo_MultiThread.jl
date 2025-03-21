
#= 
This module provides functions for MonteCarlo Integration of S and T Matrices
=#

"""
    STMonteCarloAxi_MultiThread!(SAtotal3,SAtotal4,TAtotal,SAtally3,SAtally4,TAtally,p3Max,u3MinMax,p4Max,u4MinMax,sigma,dsigmadt,Parameters,numTiterPerThread,numSiterPerThread)

# Arguments
- `SAtotal3::Array{Float64,6}` : Array of stored integration totals for S matrix for 12->34 interaction
- `SAtotal4::Array{Float64,6}` : Array of stored integration totals for S matrix for 12->43 interaction
- `TAtotal::Array{Float64,4}` : Array of stored integration totals for T matrix
- `SAtally3::Array{UInt32,5}` : Array of stored integration tallies for S matrix for 12->34 interaction
- `SAtally4::Array{UInt32,5}` : Array of stored integration tallies for S matrix for 12->43 interaction
- `TAtally::Array{UInt32,4}` : Array of stored integration tallies for T matrix
- `p3Max::Array{Float64,5}` : Array of maximum momentum values for species 3
- `u3MinMax::Array{Float64,6}` : Array of minimum and maximum theta values for species 3
- `p4Max::Array{Float64,5}` : Array of maximum momentum values for species 4
- `u4MinMax::Array{Float64,6}` : Array of minimum and maximum theta values for species 4
- `sigma::Function` : Cross section function for the interaction
- `dsigmadt::Function` : Differential cross section function for the interaction
- `Parameters::Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64,Int64,Int64}` : Tuple of parameters for the interaction
- `numTiterPerThread::Int64` : Number of T iterations per thread
- `numSiterPerThread::Int64` : Number of S iterations per thread

# Output:
- Argument arrays SAtotal,TAtotal,SAtally,TAtally are mutated to include the results of the Monte Carlo Integration.

# Calculation In Brief
- Set up worker threads
- Random Sample points in each of these domains
    - RandomPointSphere for theta and phi (for species 1,2,3,4)
    - RandomPointMomentum for p ( species 1,2 only)
- Take random points (u3,h3,p1,p2,u1,u2,h1,h2) and calculate valid p3 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Take random points (u4,h3,p1,p2,u1,u2,h1,h2) and calculate valid p4 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Update global S and T arrays with locks to prevent data races
"""
function STMonteCarloAxi_MultiThread!(SAtotal3::Array{Float64,6},SAtotal4::Array{Float64,6},TAtotal::Array{Float64,4},SAtally3::Array{UInt32,5},SAtally4::Array{UInt32,5},TAtally::Array{UInt32,4},ArrayOfLocks,p3Max::Array{Float64,5},p4Max::Array{Float64,5},u3MinMax::Array{Float64,6},u4MinMax::Array{Float64,6},sigma::Function,dsigmadt::Function,Parameters::Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64},numTiterPerThread::Int64,numSiterPerThread::Int64,p)

    # Set Parameters
    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid_st,p1_num,u1_grid_st,u1_num,p2_low,p2_up,p2_grid_st,p2_num,u2_grid_st,u2_num,p3_low,p3_up,p3_grid_st,p3_num,u3_grid_st,u3_num,p4_low,p4_up,p4_grid_st,p4_num,u4_grid_st,u4_num) = Parameters

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
    p3loc::Int64 = 0
    u3loc::Int64 = 0
    p3ploc::Int64 = 0
    u3ploc::Int64 = 0
    p4loc::Int64 = 0
    u4loc::Int64 = 0
    p4ploc::Int64 = 0
    u4ploc::Int64 = 0

    p1_grid = Grid_String_to_Type(p1_grid_st)
    p2_grid = Grid_String_to_Type(p2_grid_st)
    p3_grid = Grid_String_to_Type(p3_grid_st)
    p4_grid = Grid_String_to_Type(p4_grid_st)
    u1_grid = Grid_String_to_Type(u1_grid_st)
    u2_grid = Grid_String_to_Type(u2_grid_st)
    u3_grid = Grid_String_to_Type(u3_grid_st)
    u4_grid = Grid_String_to_Type(u4_grid_st)

    localSAtotal3 = zeros(Float64,size(SAtotal3)[1:2])
    localSAtally3 = zeros(UInt32,size(SAtally3)[1])
    localSAtotal4 = zeros(Float64,size(SAtotal4)[1:2])
    localSAtally4 = zeros(UInt32,size(SAtally4)[1])
    localp3Max = zeros(Float64,size(p3Max)[1])
    localu3Min = zeros(Float64,size(u3MinMax)[2])
    localu3Max = zeros(Float64,size(u3MinMax)[2])
    localp4Max = zeros(Float64,size(p4Max)[1])
    localu4Min = zeros(Float64,size(u4MinMax)[2])
    localu4Max = zeros(Float64,size(u4MinMax)[2])

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
        loc12 = CartesianIndex(p1loc,u1loc,p2loc,u2loc)

        fill!(localSAtally3,UInt32(0))
        fill!(localSAtally4,UInt32(0))

        if Tval != 0e0 # i.e. it is a valid interaction state

            fill!(localSAtotal3,0e0)
            fill!(localp3Max,Float64(0))
            fill!(localu3Min,Float64(0))
            fill!(localu3Max,Float64(0))
            fill!(localSAtotal4,0e0)
            fill!(localp4Max,Float64(0))
            fill!(localu4Min,Float64(0))
            fill!(localu4Max,Float64(0))
                    
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
                # For each u3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3. Therefore by sampling one u3 we are actually sampling u3 and -u3 with one or both having valid p3 states.
                #if NumStates != 0
                    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                    u3locMirror = location(u_low,u_up,u3_num,-p3v[2],u3_grid)
                    localSAtally3[u3loc] += UInt32(1)
                    localSAtally3[u3locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3loc,u3loc] += Sval
                        localp3Max[u3loc] = max(localp3Max[u3loc],p3v[1])
                        localu3Min[p3loc] = min(localu3Min[p3loc],p3v[2])
                        localu3Max[p3loc] = max(localu3Max[p3loc],p3v[2])
                    end
                end

                if NumStates == 2
                    u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3loc,u3loc] += Sval
                        localp3Max[u3loc] = max(localp3Max[u3loc],p3v[1])
                        localu3Min[p3loc] = min(localu3Min[p3loc],p3v[2])
                        localu3Max[p3loc] = max(localu3Max[p3loc],p3v[2])
                    end
                    if p3p_physical
                        p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                        Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal3[p3ploc,u3ploc] += Svalp
                        localp3Max[u3ploc] = max(localp3Max[u3ploc],p3pv[1])
                        localu3Min[p3ploc] = min(localu3Min[p3ploc],p3pv[2])
                        localu3Max[p3ploc] = max(localu3Max[p3ploc],p3pv[2])
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
                # For each u3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3. Therefore by sampling one u3 we are actually sampling u3 and -u3 with one or both having valid p3 states.
                #if NumStates != 0
                    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                    u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
                    localSAtally4[u4loc] += UInt32(1)
                    localSAtally4[u4locMirror] += UInt32(1)
                #end

                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4loc,u4loc] += Sval
                        localp4Max[u4loc] = max(localp4Max[u4loc],p4v[1])
                        localu4Min[p4loc] = min(localu4Min[p4loc],p4v[2])
                        localu4Max[p4loc] = max(localu4Max[p4loc],p4v[2])
                    end
                end

                if NumStates == 2
                    u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4loc,u4loc] += Sval
                        localp4Max[u4loc] = max(localp4Max[u4loc],p4v[1])
                        localu4Min[p4loc] = min(localu4Min[p4loc],p4v[2])
                        localu4Max[p4loc] = max(localu4Max[p4loc],p4v[2])
                    end
                    if p4p_physical
                        p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                        Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        localSAtotal4[p4ploc,u4ploc] += Svalp
                        localp4Max[u4ploc] = max(localp4Max[u4ploc],p4pv[1])
                        localu4Min[p4ploc] = min(localu4Min[p4ploc],p4pv[2])
                        localu4Max[p4ploc] = max(localu4Max[p4ploc],p4pv[2])
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
            @view(SAtally3[:,loc12]) .+= localSAtally3
            @view(SAtally4[:,loc12]) .+= localSAtally4
            if Tval != 0e0
                @view(SAtotal3[:,:,loc12]) .+= localSAtotal3
                @view(SAtotal4[:,:,loc12]) .+= localSAtotal4
                @view(p3Max[:,loc12]) .= max.(@view(p3Max[:,loc12]),localp3Max)
                @view(u3MinMax[1,:,loc12]) .= min.(@view(u3MinMax[1,:,loc12]),localu3Min)
                @view(u3MinMax[2,:,loc12]) .= max.(@view(u3MinMax[2,:,loc12]),localu3Max)
                @view(p4Max[:,loc12]) .= max.(@view(p4Max[:,loc12]),localp4Max)
                @view(u4MinMax[1,:,loc12]) .= min.(@view(u4MinMax[1,:,loc12]),localu4Min)
                @view(u4MinMax[2,:,loc12]) .= max.(@view(u4MinMax[2,:,loc12]),localu4Max)
            end 
        end 

    end # T loop

    next!(p)

    end # Thread spawn 

end # function 
