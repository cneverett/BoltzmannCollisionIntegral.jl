
#= 
This module provides functions for MonteCarlo Integration of S and T Matrices
=#

"""
    STMonteCarloAxi_Serial!(Arrays,sigma,dsigmadt,UserInputSerial)

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
- `u3MinMax::Array{Float64,6}` : Array of minimum and maximum theta values for species 4
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
- Take random points (u3,h3,p1,pu2_grid2,u1,u2,h1,h2) and calculate valid p3 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
- Take random points (u4,h3,p1,p2,u1,u2,h1,h2) and calculate valid p4 point/points 
- Find position in local S and T arrays and allocated tallies and totals accordingly.
"""
function STMonteCarloAxi_Serial!(Arrays::ScatteringArrays,sigma::Function,dsigmadt::Function,userInputSerial::BinaryUserInput)

    # Set Parameters
    Parameters = userInputSerial.Parameters
    numTiter = userInputSerial.numTiter
    numSiter = userInputSerial.numSiter
    MinMax = userInputSerial.MinMax

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

    phi1_grid_st = Parameters.phi1_grid
    phi1_num = Parameters.phi1_num

    p2_low = Parameters.p2_low
    p2_up = Parameters.p2_up
    p2_grid_st = Parameters.p2_grid
    p2_num = Parameters.p2_num

    u2_grid_st = Parameters.u2_grid
    u2_num = Parameters.u2_num

    phi2_grid_st = Parameters.phi2_grid
    phi2_num = Parameters.phi2_num

    p3_low = Parameters.p3_low
    p3_up = Parameters.p3_up
    p3_grid_st = Parameters.p3_grid
    p3_num = Parameters.p3_num

    u3_grid_st = Parameters.u3_grid
    u3_num = Parameters.u3_num

    phi3_grid_st = Parameters.phi3_grid
    phi3_num = Parameters.phi3_num

    p4_low = Parameters.p4_low
    p4_up = Parameters.p4_up
    p4_grid_st = Parameters.p4_grid
    p4_num = Parameters.p4_num

    u4_grid_st = Parameters.u4_grid
    u4_num = Parameters.u4_num

    phi4_grid_st = Parameters.phi4_grid
    phi4_num = Parameters.phi4_num

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

    p1loc::Int64 = 0
    p2loc::Int64 = 0
    u1loc::Int64 = 0
    u2loc::Int64 = 0
    h1loc::Int64 = 0
    h2loc::Int64 = 0
    p3loc::Int64 = 0
    u3loc::Int64 = 0
    h3loc::Int64 = 0
    p3ploc::Int64 = 0
    u3ploc::Int64 = 0
    h3ploc::Int64 = 0
    p4loc::Int64 = 0
    u4loc::Int64 = 0
    h4loc::Int64 = 0
    p4ploc::Int64 = 0
    u4ploc::Int64 = 0
    h4ploc::Int64 = 0
    #loc12::CartesianIndex{4} = CartesianIndex(0,0,0,0)

    p1_grid = Grid_String_to_Type(p1_grid_st)
    p2_grid = Grid_String_to_Type(p2_grid_st)
    p3_grid = Grid_String_to_Type(p3_grid_st)
    p4_grid = Grid_String_to_Type(p4_grid_st)
    u1_grid = Grid_String_to_Type(u1_grid_st)
    u2_grid = Grid_String_to_Type(u2_grid_st)
    u3_grid = Grid_String_to_Type(u3_grid_st)
    u4_grid = Grid_String_to_Type(u4_grid_st)
    h1_grid = Grid_String_to_Type(phi1_grid_st)
    h2_grid = Grid_String_to_Type(phi2_grid_st)
    h3_grid = Grid_String_to_Type(phi3_grid_st)
    h4_grid = Grid_String_to_Type(phi4_grid_st)

    p = Progress(numTiter)
    
    for _ in 1:numTiter

        # generate p1 and p2 vectors initially as to not have to re-calculate
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
        h1loc = location(phi_low,phi_up,phi1_num,p1v[3],h1_grid)
        h2loc = location(phi_low,phi_up,phi2_num,p2v[3],h2_grid)
        loc12 = CartesianIndex(p1loc,u1loc,h1loc,p2loc,u2loc,h2loc)

        SAtotalView3 = @view SAtotal3[:,:,:,loc12]
        SAtotalView4 = @view SAtotal4[:,:,:,loc12]
        SAtallyView3 = @view SAtally3[:,:,loc12]
        SAtallyView4 = @view SAtally4[:,:,loc12]

        if MinMax
            p3MaxView = @view p3Max[:,:,loc12]
            u3MinView = @view u3MinMax[1,:,:,loc12]
            u3MaxView = @view u3MinMax[2,:,:,loc12]
            p4MaxView = @view p4Max[:,:,loc12]
            u4MinView = @view u4MinMax[1,:,:,loc12]
            u4MaxView = @view u4MinMax[2,:,:,loc12]
        end
        
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
                # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). Therefore by sampling one u3,h3 we are actually sampling u3 and -u3 and h3, mod(h3+1,2) with one or both having valid p3 states.
                #if NumStates != 0
                    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
                    u3locMirror = location(u_low,u_up,u3_num,-p3v[2],u3_grid)
                    h3loc = location(phi_low,phi_up,phi3_num,p3v[3],h3_grid)
                    h3locMirror = location(phi_low,phi_up,phi3_num,mod(p3v[3]+1e0,2e0),h3_grid)
                    SAtallyView3[u3loc,h3loc] += UInt32(1)
                    SAtallyView3[u3locMirror,h3locMirror] += UInt32(1)
                #end
  
                # Calculate S Array totals
                if NumStates == 1
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView3[p3loc,u3loc,h3loc] += Sval
                        if MinMax
                            p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                            u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                            u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
                        end
                    end
                end

                if NumStates == 2
                    if p3_physical
                        p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                        Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView3[p3loc,u3loc,h3loc] += Sval
                        if MinMax
                            p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                            u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                            u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
                        end
                    end
                    if p3p_physical
                        u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                        h3ploc = location(phi_low,phi_up,phi3_num,p3pv[3],h3_grid)
                        p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                        Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView3[p3ploc,u3ploc,h3ploc] += Svalp
                        if MinMax
                            p3MaxView[u3ploc,h3ploc] = max(p3MaxView[u3ploc,h3ploc],p3pv[1])
                            u3MinView[p3ploc,h3ploc] = min(u3MinView[p3ploc,h3ploc],p3pv[2])
                            u3MaxView[p3plo,h3ploc] = max(u3MaxView[p3ploc,h3ploc],p3pv[2])
                        end
                    end
                end

            # === p4 === #
                #set random p4 direction 
                p4v .= pv
                p4pv .= pv

                # Calculate p4 value
                (p4_physical,p4p_physical,NumStates) = Momentum3Value!(p4v,p4pv,p2v,p1v,mu2,mu1,mu4,mu3)

                # S Array Tallies
                # For each u3,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u3 and a shift in h4 by pi i.e. Mod(h4+1,2). Therefore by sampling one u3 we are actually sampling u3/h4 and -u3/mod(h4+1,2) with one or both having valid p4 states.
                #if NumStates != 0
                    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
                    h4loc = location(phi_low,phi_up,phi4_num,p4v[3],h4_grid)
                    u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
                    h4locMirror = location(phi_low,phi_up,phi4_num,mod(p4v[3]+1e0,2e0),h4_grid)
                    SAtallyView4[u4loc,h4loc] += UInt32(1)
                    SAtallyView4[u4locMirror,h4locMirror] += UInt32(1)
                #end
    
                # Calculate S Array totals
                if NumStates == 1
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView4[p4loc,u4loc,h4loc] += Sval
                        if MinMax
                            p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                            u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                            u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
                        end
                    end
                end

                if NumStates == 2
                    if p4_physical
                        p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                        Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView4[p4loc,u4loc,h4loc] += Sval
                        if MinMax
                            p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                            u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                            u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
                        end
                    end
                    if p4p_physical
                        u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                        h4ploc = location(phi_low,phi_up,phi4_num,p4pv[3],h4_grid)
                        p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                        Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                        SAtotalView4[p4ploc,u4ploc,h4ploc] += Svalp
                        if MinMax
                            p4MaxView[u4ploc,h4ploc] = max(p4MaxView[u4ploc,h4ploc],p4pv[1])
                            u4MinView[p4ploc,h4ploc] = min(u4MinView[p4ploc,h4ploc],p4pv[2])
                            u4MaxView[p4ploc,h4ploc] = max(u4MaxView[p4ploc,h4ploc],p4pv[2])
                        end
                    end
                end

            end # Sloop =#

        else # no valid interaction state
            # add one to tally of all relevant S tallies i.e. all momenta and all angles as no emission states are possible
            SAtallyView3 .+= UInt32(1)
            SAtallyView4 .+= UInt32(1)
        end

        # assign to T arrays
        TAtotal[loc12] += Tval # ST[3] doesn't change with S loop
        TAtally[loc12] += UInt32(1)

        next!(p)

    end # Tloop

    finish!(p)

    return nothing

end