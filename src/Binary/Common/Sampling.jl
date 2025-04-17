"""
    UniformSampling()

Uniform MC sampling of outgoing momentum states in the un-boosted Lab frame.
"""
function UniformSampling!(pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p3pv::Vector{Float64},p4pv::Vector{Float64},#=Sval::Float64,Svalp::Float64,=#dsigmadt::Function,SAtotalView3::Array{Float64,2},SAtotalView4::Array{Float64,2},SAtallyView3::Array{Float64,3},SAtallyView4::Array{Float64,3},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # Unpack locations, NOT NEEDED
    #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate random p direction in Lab frame for use in both p3 and p4 calculations
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
        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
        h3locMirror = location(h_low,h_up,h3_num,mod(p3v[3]+1e0,2e0),h3_grid)
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
            h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
            p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
            Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView3[p3ploc,u3ploc,h3ploc] += Svalp
            if MinMax
                p3MaxView[u3ploc,h3ploc] = max(p3MaxView[u3ploc,h3ploc],p3pv[1])
                u3MinView[p3ploc,h3ploc] = min(u3MinView[p3ploc,h3ploc],p3pv[2])
                u3MaxView[p3ploc,h3ploc] = max(u3MaxView[p3ploc,h3ploc],p3pv[2])
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
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
    h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
    SAtallyView4[u4loc,h4loc] += UInt32(1)
    SAtallyView4[u4locMirror,h4locMirror] += UInt32(1)


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
            h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
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

end

function UniformSampling2!(pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p3pv::Vector{Float64},p4pv::Vector{Float64},#=Sval::Float64,Svalp::Float64,=#dsigmadt::Function,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # Unpack locations, NOT NEEDED
    #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate random p direction in Lab frame for use in both p3 and p4 calculations
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
        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
        h3locMirror = location(h_low,h_up,h3_num,mod(p3v[3]+1e0,2e0),h3_grid)
        SAtallyView3[u3loc,h3loc] += UInt32(1)
        SAtallyView3[u3locMirror,h3locMirror] += UInt32(1)
    #end

    # Calculate S Array totals
    if NumStates == 1
        if p3_physical
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView3[p3loc,u3loc,h3loc] += Sval

            pVector!(p4v,p3v,p1v,p2v)
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
            h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
            Svalp = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtallyView4[u4loc,h4loc] += UInt32(1)
            SAtotalView4[p4loc,u4loc,h4loc] += Svalp


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

            pVector!(p4v,p3v,p1v,p2v)
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
            h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
            Svalp = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtallyView4[u4loc,h4loc] += UInt32(1)
            SAtotalView4[p4loc,u4loc,h4loc] += Svalp

            if MinMax
                p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
            end
        end
        if p3p_physical
            u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
            h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
            p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
            Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView3[p3ploc,u3ploc,h3ploc] += Svalp

            pVector!(p4v,p3pv,p1v,p2v)
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
            h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
            Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtallyView4[u4loc,h4loc] += UInt32(1)
            SAtotalView4[p4loc,u4loc,h4loc] += Sval

            if MinMax
                p3MaxView[u3ploc,h3ploc] = max(p3MaxView[u3ploc,h3ploc],p3pv[1])
                u3MinView[p3ploc,h3ploc] = min(u3MinView[p3ploc,h3ploc],p3pv[2])
                u3MaxView[p3ploc,h3ploc] = max(u3MaxView[p3ploc,h3ploc],p3pv[2])
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
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    u4locMirror = location(u_low,u_up,u4_num,-p4v[2],u4_grid)
    h4locMirror = location(h_low,h_up,h4_num,mod(p4v[3]+1e0,2e0),h4_grid)
    SAtallyView4[u4loc,h4loc] += UInt32(1)
    SAtallyView4[u4locMirror,h4locMirror] += UInt32(1)


    # Calculate S Array totals
    if NumStates == 1
        if p4_physical
            p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
            Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView4[p4loc,u4loc,h4loc] += Sval

            pVector!(p3v,p4v,p1v,p2v)
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
            h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
            Svalp = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtallyView3[u3loc,h3loc] += UInt32(1)
            SAtotalView3[p3loc,u3loc,h3loc] += Svalp

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

            pVector!(p3v,p4v,p1v,p2v)
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
            h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
            Svalp = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtallyView3[u3loc,h3loc] += UInt32(1)
            SAtotalView3[p3loc,u3loc,h3loc] += Svalp
 
            if MinMax
                p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
            end
        end
        if p4p_physical
            u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
            h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
            p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
            Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtotalView4[p4ploc,u4ploc,h4ploc] += Svalp

            pVector!(p3v,p4v,p1v,p2v)
            p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
            u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
            h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
            Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
            SAtallyView3[u3loc,h3loc] += UInt32(1)
            SAtotalView3[p3loc,u3loc,h3loc] += Sval

            if MinMax
                p4MaxView[u4ploc,h4ploc] = max(p4MaxView[u4ploc,h4ploc],p4pv[1])
                u4MinView[p4ploc,h4ploc] = min(u4MinView[p4ploc,h4ploc],p4pv[2])
                u4MaxView[p4ploc,h4ploc] = max(u4MaxView[p4ploc,h4ploc],p4pv[2])
            end
        end
    end

end

"""
    ImportanceSampling()

Importance MC sampling of outgoing momentum states in the boosted centre-of-momentum frame. i.e. importance sampling weighted by the de-boosting of angles from the COM frame to the Lab frame.
"""
function ImportanceSampling(pv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p3pv::Vector{Float64},p4pv::Vector{Float64},βv::Vector{Float64},pCv::Vector{Float64}#=,Sval::Float64,Svalp::Float64=#,dsigmadt::Function,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

        # Unpack parameters
        (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

        # Unpack locations NOT NEEDED
        #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate COM frame velocity vector 
    betaVec!(βv,p1v,p2v,mu1,mu2)
    # generate random p direction in COM frame and boost back to Lab frame for use in both p3 and p4 calculations
    RPointSphereCosThetaPhi!(pCv)
    RPointSphereBoost!(pv,βv,pCv)

    # === p3 === #
    #set random p3 direction 
    p3v .= pv
    p3pv .= pv

    # Calculate p3 value
    (p3_physical,p3p_physical,NumStates) = Momentum3Value!(p3v,p3pv,p1v,p2v,mu1,mu2,mu3,mu4)

    # S Array Tallies
    # For each u3,h3 sampled, p3 will be + or -ve, corresponding to a change in sign of u3 and a rotation of h3 by pi i.e. mod(h3+1,2). As this sampling is conducted in the COM frame, it would give the incorrect statistics to include the mirrored states. Therefore only the un-mirrored states are included.
    # un-mirrored states have sign(p3v[2]) == sign(pv[2]) 

    #if NumStates != 0
        u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
        h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
        SAtallyView3[u3loc,h3loc] += UInt32(1)
    #end

    # Calculate S Array totals
    if NumStates == 1
        if p3_physical
            if sign(p3v[2]) == sign(pv[2])
                p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView3[p3loc,u3loc,h3loc] += Sval/prob

                if Sval == 0e0
                    error("Sval = 0e0")
                end

                if MinMax
                    p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                    u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                    u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
                end
            end
        end
    end

    if NumStates == 2
        if p3_physical
            if sign(p3v[2]) == sign(pv[2])
                p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
                Sval = SValue3(p3v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView3[p3loc,u3loc,h3loc] += Sval/prob

                if Sval == 0e0
                    error("Sval = 0e0")
                end

                if MinMax
                    p3MaxView[u3loc,h3loc] = max(p3MaxView[u3loc,h3loc],p3v[1])
                    u3MinView[p3loc,h3loc] = min(u3MinView[p3loc,h3loc],p3v[2])
                    u3MaxView[p3loc,h3loc] = max(u3MaxView[p3loc,h3loc],p3v[2])
                end
            end
        end
        if p3p_physical
            if sign(p3pv[2]) == sign(pv[2])
                u3ploc = location(u_low,u_up,u3_num,p3pv[2],u3_grid)
                h3ploc = location(h_low,h_up,h3_num,p3pv[3],h3_grid)
                p3ploc = location(p3_low,p3_up,p3_num,p3pv[1],p3_grid)
                Svalp = SValue3(p3pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView3[p3ploc,u3ploc,h3ploc] += Svalp/prob

                if Sval == 0e0
                    error("Sval = 0e0")
                end

                if MinMax
                    p3MaxView[u3ploc,h3ploc] = max(p3MaxView[u3ploc,h3ploc],p3pv[1])
                    u3MinView[p3ploc,h3ploc] = min(u3MinView[p3ploc,h3ploc],p3pv[2])
                    u3MaxView[p3ploc,h3ploc] = max(u3MaxView[p3ploc,h3ploc],p3pv[2])
                end
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
    # For each u4,h4 sampled, p4 will be + or -ve, corresponding to a change in sign of u4 and a rotation of h4 by pi i.e. mod(h4+1,2). As this sampling is conducted in the COM frame, it would give the incorrect statistics to include the mirrored states. Therefore only the un-mirrored states are included.
    # un-mirrored states have sign(p4v[2]) == sign(pv[2]) 
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    SAtallyView4[u4loc,h4loc] += UInt32(1)

    # Calculate S Array totals
    if NumStates == 1
        if p4_physical
            if sign(p4v[2]) == sign(pv[2])
                p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView4[p4loc,u4loc,h4loc] += Sval/prob
                if MinMax
                    p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                    u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                    u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
                end
            end
        end
    end

    if NumStates == 2
        if p4_physical
            if sign(p4v[2]) == sign(pv[2])
                p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
                Sval = SValue4(p4v,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView4[p4loc,u4loc,h4loc] += Sval/prob
                if MinMax
                    p4MaxView[u4loc,h4loc] = max(p4MaxView[u4loc,h4loc],p4v[1])
                    u4MinView[p4loc,h4loc] = min(u4MinView[p4loc,h4loc],p4v[2])
                    u4MaxView[p4loc,h4loc] = max(u4MaxView[p4loc,h4loc],p4v[2])
                end
            end
        end
        if p4p_physical
            if sign(p4pv[2]) == sign(pv[2])
                u4ploc = location(u_low,u_up,u4_num,p4pv[2],u4_grid)
                h4ploc = location(h_low,h_up,h4_num,p4pv[3],h4_grid)
                p4ploc = location(p4_low,p4_up,p4_num,p4pv[1],p4_grid)
                Svalp = SValue4(p4pv,p1v,p2v,dsigmadt,mu1,mu2,mu3,mu4)
                prob = pdfBoost(βv[1],pCv[1])
                SAtotalView4[p4ploc,u4ploc,h4ploc] += Svalp/prob
                if MinMax
                    p4MaxView[u4ploc,h4ploc] = max(p4MaxView[u4ploc,h4ploc],p4pv[1])
                    u4MinView[p4ploc,h4ploc] = min(u4MinView[p4ploc,h4ploc],p4pv[2])
                    u4MaxView[p4ploc,h4ploc] = max(u4MaxView[p4ploc,h4ploc],p4pv[2])
                end
            end
        end
    end

    #if (sum(SAtotalView3) == 0e0) || (sum(SAtotalView4) == 0e0)
    #    println("No valid p3 or p4 states found")
    #end

    return nothing


end

function ImportanceSampling2(p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p1cv::Vector{Float64},p1cBv::Vector{Float64},Tval::Float64,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # Unpack locations NOT NEEDED
    #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate importance sampled points 
    importance!(p3v,p4v,p1v,p2v,p1cv,p1cBv,mu1,mu2,mu3,mu4)
    prob3 = p3v[5]
    prob4 = p4v[5]

    p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
    p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    SAtallyView3[u3loc,h3loc] += UInt32(1)
    SAtallyView4[u4loc,h4loc] += UInt32(1)

    Sval3,Sval4 = SValueTest(p3v,p4v,Tval,p1v,p2v,mu1,mu2,mu3,mu4)
    SAtotalView3[p3loc,u3loc,h3loc] += Sval3*prob3
    SAtotalView4[p4loc,u4loc,h4loc] += Sval4*prob4

    #=if p4v[1] <= 1e-8
        println("p4 = $(p4v[1])")
        println("sp = $(Sval4*prob4)")
        println("p4loc = $p4loc")

    end=#

    #println("prob = $prob3")
    #println("sval3 = $Sval3")

    #if (sum(SAtotalView3) == 0e0) || (sum(SAtotalView4) == 0e0)
    #    println("No valid p3 or p4 states found")
    #end

    return nothing


end

function ImportanceSampling3(p1v::Vector{Float64},p2v::Vector{Float64},p3v::Vector{Float64},p4v::Vector{Float64},p1cv::Vector{Float64},p1cBv::Vector{Float64},dsigmadt::Function,SAtotalView3::AbstractArray{Float64,3},SAtotalView4::AbstractArray{Float64,3},SAtallyView3::AbstractArray{UInt32,2},SAtallyView4::AbstractArray{UInt32,2},Parameters::Tuple{Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Int64,GridType,Int64,GridType,Int64,GridType, Float64,Float64,Float64,Float64},#=Locations::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64},=#MinMax::Bool;p3MaxView=nothing,u3MinView=nothing,u3MaxView=nothing,p4MaxView=nothing,u4MinView=nothing,u4MaxView=nothing)

    # Unpack parameters
    (p3_low,p3_up,p3_num,p3_grid,u3_num,u3_grid,h3_num,h3_grid,p4_low,p4_up,p4_num,p4_grid,u4_num,u4_grid,h4_num,h4_grid,mu1,mu2,mu3,mu4) = Parameters

    # Unpack locations NOT NEEDED
    #(p3loc,u3loc,u3locMirror,h3loc,h3locMirror,p3ploc,u3ploc,h3ploc,p4loc,u4loc,u4locMirror,h4loc,h4locMirror,p4ploc,u4ploc,h4ploc) = Locations

    # generate importance sampled points 
    tSmol,sSmol = importance!(p3v,p4v,p1v,p2v,p1cv,p1cBv,mu1,mu2,mu3,mu4)
    prob3 = p3v[5]
    prob4 = p4v[5]

    p3loc = location(p3_low,p3_up,p3_num,p3v[1],p3_grid)
    p4loc = location(p4_low,p4_up,p4_num,p4v[1],p4_grid)
    u3loc = location(u_low,u_up,u3_num,p3v[2],u3_grid)
    h3loc = location(h_low,h_up,h3_num,p3v[3],h3_grid)
    u4loc = location(u_low,u_up,u4_num,p4v[2],u4_grid)
    h4loc = location(h_low,h_up,h4_num,p4v[3],h4_grid)
    SAtallyView3[u3loc,h3loc] += UInt32(1)
    SAtallyView4[u4loc,h4loc] += UInt32(1)

    Sval3,Sval4 = SValueTest2(p3v,p4v,dsigmadt,tSmol,sSmol,p1v,p2v,mu1,mu2,mu3,mu4)
    SAtotalView3[p3loc,u3loc,h3loc] += Sval3*prob3
    SAtotalView4[p4loc,u4loc,h4loc] += Sval4*prob4

    #=if p4v[1] <= 1e-8
        println("p4 = $(p4v[1])")
        println("sp = $(Sval4*prob4)")
        println("p4loc = $p4loc")

    end=#

    #println("prob = $prob3")
    #println("sval3 = $Sval3")

    #if (sum(SAtotalView3) == 0e0) || (sum(SAtotalView4) == 0e0)
    #    println("No valid p3 or p4 states found")
    #end

    return nothing


end


function Inv_cdfElePhoElePho(sS::Float64)

    # generates a random number x between [0,1] and returns the inverse cdf of the Compton cross section i.e. tSmol value.

    #cdf(tSmol) is
    #(((tS + sS*(sS + tS)) * (8*tS + sS*(16 + sS^3 + 2*sS^2*(9 + tS) + tS*(16 + tS) + sS*(32 + tS*(11 + tS)))))/(2*(1 + sS)^2 *(sS + tS)) + sS*(-8 + (-4 + sS)*sS)*log(1 + sS + tS + tS/sS))/(sS^3*(8/sS + (sS*(2 + sS))/(2*(1 + sS)^2) + ((-8 + (-4 + sS)*sS)*log(1 + sS))/sS^2))

    # tSmol can be normalised by sS^2/(sS+1) to be tSmolN, then cos(ts) = 1-2tSmolN.
    # Therefore prob(cos(ts)) = |-2*prob(tSmolN)|

    x = rand(Float64)

    #logsS=log10(sS)
    logsS=log10(sS)

    # Inverse cdf(tSmol) can be approximate by to regions
    if logsS < 1e0 #0.8

        tSmolN = (x - 1)
        #tS *= sS^2/(sS + 1)

        # more floating point safe than acos()
        ts = 2*atan(sqrt(-tSmolN),sqrt(tSmolN+1))/pi
        
        prob_cts = 1e0 #2e0 # factor of 2 absorbed into cos(theta) probability distribution

    else #logsS >= log10

        tSmolN = x^(logsS)-1

        # more floating point safe than acos()
        ts = 2*atan(sqrt(-tSmolN),sqrt(tSmolN+1))/pi

        prob = (1+tSmolN)^(1/logsS-1) 
        prob /= logsS 

        prob_cts = prob #2e0 # factor of 2 absorbed into cos(theta) probability distribution

    end

    tSmol = tSmolN * sS^2/(sS + 1)

    if tS < -1.0 #tS < -sS^2/(sS+1)
        println("tS = $tS")
        println("sS = $sS")
        println("tS min = $(-sS^2/(sS+1))")
        error("tS < sS^2/(sS+1)")
    end

    #println("tS = $(tS*sS^2/(sS+1))")

    return (tSmol,ts,prob_cts)

end

function RPointPhi()

    phi = 2*rand(Float64)

    return phi

end

function p1cv!(p1cv::Vector{Float64},p1cBv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},m1::Float64,m2::Float64,sS::Float64)

    p1 = p1v[1]
    ct1 = p1v[2]
    st1 = sqrt(1e0 - ct1^2)

    p2 = p2v[1]
    ct2 = p2v[2]
    st2 = sqrt(1e0 - ct2^2)

    sh1h2, ch1h2 = sincospi(p1v[3]-p2v[3])

    m12 = m2^2
    m22 = m2^2

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    E1::Float64 = Es1 + m1
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    E2::Float64 = Es2 + m2

    val1 = (ct1*ct2+ch1h2*st1*st2)
    val2 = 2*p1*p2*val1
    val3 = val2 + p1^2 + p2^2

    w = asinh(sqrt(val3)/sqrt(sS+(m1+m2)^2))

    # calculate p1cv value
    if (y=p2/p1) <= 1e-6
        t1c = 2*atan((1-val1^2)*y^2/8)/pi
    elseif y >= 1e6
        t1c = 2*atan((1-val1)/(2*(1+val1)) + (val1-1)/(1+val1)/y)/pi
    else
        y = p1+p2*val1
        y /= sqrt(val3)
        t1c = acos(y)/pi
    end

    p1cv[2] = cospi(t1c)
    p1cv[4] = t1c

    y = sh1h2*st1*st2*sqrt(val3)
    x = ct2*(p1+p2*val1)-ct1*(p2+p1*val1)
    p1cv[3] = mod(atan(y,x)/pi,2)

    # Boost p1c to CM frame
    LorentzBoost!(p1cBv,p1cv,m1,w)

    return nothing

end

function tstar(tS::Float64,sS::Float64,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # calculate tstar value

    #cts = (sS+(m1+m2)^2)*(2*tS+sS+(2*m12+m32-m42-4*m1*m3+2*m1*m2)) + (m12-m22)*(m32-m42)
    #cts /= 4*F12*F34

    #cts= ((m1 - m2)*(m1 + m2)*(m3 - m4)*(m3 + m4) + ((m1 + m2)^2 + sS)*(2*m1*(m1 + m2) - 4*m1*m3 + m3^2 - m4^2 + sS + 2*tS))/(sqrt(sS*(4*m1* m2 + sS))*sqrt((m3^2 - m4^2)^2 - 2*(m3^2 + m4^2)*((m1 + m2)^2 + sS) + ((m1 + m2)^2 + sS)^2))

    # cts with m values put in for compton
        #cts = 2*tS/sS^2 + 2*tS/sS + 1 
    # this is potentially very small negative +1 so subject to float issues
    # use 1-cts = 2*sin^2(ts/2) = -2*tS*(sS+1)/sS^2
    # use 1+cts = 2*cos^2(ts/2) = 2(tS*(sS+1)+sS^2)/sS^2
        #onemcts = -2*tS*(sS+1)/sS^2
        #ts = 2*asin(sqrt(onemcts/2))   
        #ts = 2*atan(-tS*(sS+1),tS*(sS+1)+sS^2)/pi

        # tS normalised by sS^2/(sS+1)
        ts = 2*atan(sqrt(-tS),sqrt(tS+1))/pi

    #println(tS)
    #println(sS)
    #println(sS^2/(sS+1))
    #=if abs(cts) > 1e0
        println("tS=$tS")
        println("sS=$sS")
        println("tS max = $(sS^2/(sS+1))")
        println("cts = $cts")
        error("cts > 1e0")
    end=#
    #if cts > 1e0
    #    cts = 1e0
    #elseif cts < -1e0
    #    cts = -1e0
    #end

    return ts 

end

function p3p4v!(p3v::Vector{Float64},p4v::Vector{Float64},p1cBv::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},ts::Float64,hs::Float64,sS::Float64,m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    t1cB = p1cBv[4]
    (st1cB,ct1cB) = sincospi(t1cB)

    (sts,cts) = sincospi(ts)

    p1 = p1v[1]
    p2 = p2v[1]

    ct1 = p1v[2]
    ct2 = p2v[2]
    st1 = sqrt(1-ct1^2)
    st2 = sqrt(1-ct2^2)

    (sh1,ch1) = sincospi(p1v[3])
    (sh2,ch2) = sincospi(p2v[3])
    ch1h2 = cospi(p1v[3]-p2v[3])

    m12 = m1^2
    m22 = m2^2

    Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
    E1::Float64 = Es1 + m1
    Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
    E2::Float64 = Es2 + m2

    val1 = (ct1*ct2+ch1h2*st1*st2)
    val2 = 2*p1*p2*val1
    val3 = val2 + p1^2 + p2^2
    val3 = 2*p1*p2*(val1-1) + (p1+p2)^2

    w = asinh(sqrt(val3)/sqrt(sS+(m1+m2)^2))

    # calculate ct3cs and ct4cs
    #t3cs = acos(cts*ct1cs+cospi(hs)*sts*st1cs)/pi
    if ts <= 1e-5 && t1cB <= 1e-5
        t3cB = sqrt(-2*cospi(hs)*ts*t1cB+ts^2+t1cB^2)/pi
    else
        t3cB = acos(cts*ct1cB+cospi(hs)*sts*st1cB)/pi
    end

    # calculate h3,h4 from ts, hs
    A = atan(sts*sinpi(hs),cts*st1cB-cospi(hs)*sts*ct1cB)/pi
    #A1 = atan((sts*sinpi(hs))/(cts*st1cB-cospi(hs)*sts*ct1cB))/pi
    #A2 = acos((cts*st1cB-cospi(hs)*sts*ct1cB)/sinpi(t3cB))/pi
    h3cB = mod(p1cBv[3] + A,2)

    # Boosted frame angles and momenta
    p3v[4] = t3cB
    p4v[4] = 1-t3cB
    p3v[3] = h3cB
    p4v[3] = mod(h3cB+1e0,2)
    p3B = InvariantFluxSmall(sS,m3,m4)/sqrt(sS+(m1+m2)^2)
    p4B = p3B
    p3v[1] = p3B
    p4v[1] = p4B

    # Doppler Factors
    DopplerFactor!(p3v,m3,-w)
    DopplerFactor!(p4v,m4,-w)

    # De-boost p3 and p4 (modifies p3v and p4v)
    LorentzBoost!(p3v,p3v,m3,-w)
    LorentzBoost!(p4v,p4v,m4,-w)

    # doppler factor
    #p3v[4] /= (γ+γβ*cospi(t3cs))^2 
    #p3v[4] = 1/(exp(w)*cospi(t3cB/2)^2+exp(-w)*sinpi(t3cB/2)^2)^2
    #p4v[4] /= (γ-γβ*cospi(t3cs))^2
    #p4v[4] = 1/(exp(w)*sinpi(t3cB/2)^2+exp(-w)*cospi(t3cB/2)^2)^2

    # p3 value
    #m32 = m3^2
    #m42 = m4^2
    #p3s = InvariantFluxSmall(sS,m3,m4)/sqrt(sS+(m1+m2)^2)
    #E3s = sqrt(p3s^2+m32)
    #testE3s = 1/2*(sS+(m12+m32+2*m1*m2+m22-m42))/sqrt(sS+(m1+m2)^2)
    #p3s = sqrt(E3s^2-m32)
    
    #p3v[1] = γ*p3s*ct3cs+γβ*E3s
    #p3c = p3s*(exp(w)*cospi(t3cB/2)^2+exp(-w)*sinpi(t3cB/2)^2)+γβ*m32/(E3s+p3s)
    #p3c = p3c/ct3c*sign(ct3c)
    #p3v[1] = p3c 

    # p4 value
    #p4s = p3s
    #E4s = sqrt(p4s^2+m42)

    #p4v[1] = -γ*p4s*ct3cs+γβ*E4s
    #p4c = p4s*(exp(w)*sinpi(t3cs/2)^2+exp(-w)*cospi(t3cs/2)^2)+γβ*m42/(E4s+p4s)
    #p4c = p4c/ct4c#*sign(ct4c)
    #p4v[1] = p4c

    # rotate to lab frame
    βc = sqrt(val3)#/(E1+E2) 
    tβ = acos((p1*ct1 + p2*ct2)/βc)/pi #/ (E1+E2) # add an approximation here for p1/p2 << 1 etc
    if (y=p1/p2) <= 1e-6
        tβ = acos(ct2 + st2*(-ch1h2*ct2*st1 + ct1*st2)*y)/pi
    elseif y >= 1e6
        tβ = acos(ct1 + st1*(ct2*st1 - ch1h2*ct1*st2)/y)/pi
    else
        tβ = acos((p1*ct1 + p2*ct2)/βc)/pi
    end
    (stβ,ctβ) = sincospi(tβ)
    x = p1*st1*ch1 + p2*st2*ch2
    y = p1*st1*sh1 + p2*st2*sh2
    hβ = mod(atan(y,x)/pi,2)
    (shβ,chβ) = sincospi(hβ)
    (st3c,ct3c) = sincospi(p3v[4])
    (st4c,ct4c) = sincospi(p4v[4])
    (sh3c,ch3c) = sincospi(p3v[3])
    (sh4c,ch4c) = sincospi(p4v[3])

    p3v[4] = acos(ct3c*ctβ + st3c*stβ*ch3c)/pi
    p3v[2] = cospi(p3v[4])
    x = -st3c*sh3c*shβ + chβ*(st3c*ch3c*ctβ + ct3c*stβ)
    y = st3c*sh3c*chβ + shβ*(st3c*ch3c*ctβ + ct3c*stβ)
    p3v[3] = mod(atan(y,x)/pi,2)

    p4v[4] = acos(ct4c*ctβ + st4c*stβ*ch4c)/pi
    p4v[2] = cospi(p4v[4])
    x = -st4c*sh4c*shβ + chβ*(st4c*ch4c*ctβ + ct4c*stβ)
    y = st4c*sh4c*chβ + shβ*(st4c*ch4c*ctβ + ct4c*stβ)
    p4v[3] = mod(atan(y,x)/pi,2)

    # energy check
    eng_error = 1-(sqrt(m3^2+p3v[1]^2)+sqrt(m4^2+p4v[1]^2))/(E1+E2)
    mom_z_error = 1-(p3v[1]*p3v[2]+p4v[1]*p4v[2])/(p1v[1]*p1v[2]+p2v[1]*p2v[2])
    mom_x_error = 1-(p3v[1]*sin(acos(p3v[2]))*cospi(p3v[3])+p4v[1]*sin(acos(p4v[2]))*cospi(p4v[3]))/(p1v[1]*sin(acos(p1v[2]))*cospi(p1v[3])+p2v[1]*sin(acos(p2v[2]))*cospi(p2v[3]))
    #mom_x_error = 1-(p3v[1]*sinpi(p3v[4])*cospi(p3v[3])+p4v[1]*sinpi(p4v[4])*cospi(p4v[3]))/(p1v[1]*sin(acos(p1v[2]))*cospi(p1v[3])+p2v[1]*sin(acos(p2v[2]))*cospi(p2v[3]))
    if eng_error > 1e-4 || mom_z_error > 1e-4 #|| mom_x_error > 1e-2
        println("eng_error = $eng_error")
        println("mom_z_error = $mom_z_error")
        println("mom_x_error = $mom_x_error")
        println("w=$w")
        println("p1 = $(p1v[1])")
        println("ct1 = $(p1v[2])")
        println("h1 = $(p1v[3])")
        println("p2 = $(p2v[1])")
        println("ct2 = $(p2v[2])")
        println("h2 = $(p2v[3])")
        println("p3 = $(p3v[1])")
        println("ct3 = $(p3v[2])")
        println("h3 = $(p3v[3])")
        println("p4 = $(p4v[1])")
        println("ct4 = $(p4v[2])")
        println("h4 = $(p4v[3])")
        println("")
    end

    #=if p1 >= 1e4 && p2 <= 1e0 && log10(sS) > 4.0
        println("ss = $sS")
        println("t1 = $(acos(p1v[2])/pi)")
        println("t2 = $(acos(p2v[2])/pi)")
        println("t1c = $(acos(p1cv[2])/pi)")
        println("t1cs = $t1cs")
        println("ts = $ts")
        println("t3cs = $t3cs")
        println("t3 = $(acos(p3v[2])/pi)")
        println("t4 = $(acos(p4v[2])/pi)")
        println("")
    end=#


end

function importance!(p3v::Vector{Float64},p4v::Vector{Float64},p1v::Vector{Float64},p2v::Vector{Float64},p1cv::Vector{Float64},p1cBv::Vector{Float64},m1::Float64,m2::Float64,m3::Float64,m4::Float64)

    # calculate random tS and hs and then ts

        # pre-defining terms for efficiency 
        p1::Float64 = p1v[1]
        p2::Float64 = p2v[1]
    
        ct1::Float64 = p1v[2] #cospi(p1v[2])
        ct2::Float64 = p2v[2] #cospi(p2v[2]) 
    
        st1::Float64 = sqrt(1e0-p1v[2]^2) #sinpi(p1v[2])
        st2::Float64 = sqrt(1e0-p2v[2]^2) #sinpi(p2v[2])
    
        ch1h2::Float64 = cospi(p1v[3]-p2v[3])
    
        m32 = m3^2
        m42 = m4^2
        m12 = m1^2
        m22 = m2^2
        
        Es1::Float64 = m1 != 0e0 ? (p1^2)/(sqrt(m12+p1^2)+m1) : p1
        Es1s::Float64 = Es1/p1
        E1::Float64 = Es1 + m1
     
        Es2::Float64 = m2 != 0e0 ? (p2^2)/(sqrt(m22+p2^2)+m2) : p2
        Es2s::Float64 = Es2/p2
        E2::Float64 = Es2 + m2

    sBig::Float64 = (m1+m2)^2
    #sSmol::Float64 = 2*(m1*Es2 + m2*Es1 + Es1*Es2 - p1*p2*(ct1*ct2+ch1h2*st1*st2))
    sSmol::Float64 = 2*p1*p2*(-ct1*ct2 -ch1h2*st1*st2 + Es1s*Es2s + m1*Es2s/p1 + m2*Es1s/p2)

    (tSmol,ts,prob) = Inv_cdfElePhoElePho(sSmol)
    hs = RPointPhi()
    #ts = tstar(tSmol,sSmol,m1,m2,m3,m4)

    # calculate p1cv
    p1cv!(p1cv,p1cBv,p1v,p2v,m1,m2,sSmol)

    # calculate p3v and p4v
    p3p4v!(p3v,p4v,p1cv,p1v,p2v,ts,hs,sSmol,m1,m2,m3,m4)

    p3v[5] /= prob
    p4v[5] /= prob

    return (tSmol,sSmol)

end


"""
    LorentzBoost!(pBv::Vector{Float64},pv::Vector{Float64},m::Float64,w::Float64)

    Calculates the Lorentz boost (in z direction) components of the momentum vector `pv` of a particle with (normalised) mass `m` by the rapidity `w` (can be positive or negative), placing the results in the boosted vector `pBv``. 

    Momentum vectors have coponents [p,cos(theta),phi,theta]
"""
function LorentzBoost!(pBv::Vector{Float64},pv::Vector{Float64},m::Float64,w::Float64)

    p::Float64 = pv[1]
    t::Float64 = pv[4]
    h::Float64 = pv[3]

    hB::Float64 = h # phi unaffected by boost

    if m == 0e0

        ew = exp(w)

        pB = p*(ew*sinpi(t/2)^2+cospi(t/2)^2/ew)
        tB = 2*atan(ew*tanpi(t/2))/pi

    elseif p/m <= 1e-4 # small p approximation 
        # valid to at least second order in p/m

        pB = sign(w)*(m*sinh(w)-p*cospi(t)*cosh(w))
        tB = 2*atan(p*sinpi(t),2(sign(w)*sinh(w)))/pi

    elseif p/m >= 1e4 # large p approximation
        # valid to at least second order in m/p

        ew = exp(w)
        ttd2 = tanpi(t/2)

        pB = p*(ew*sinpi(t/2)^2+cospi(t/2)^2/ew)
        pB += (m^2/(2*p))*sinh(w)*(sinpi(t/2)^2*ew-cospi(t/2)^2/ew)/(sinpi(t/2)^2*ew+cospi(t/2)^2/ew)

        tB = 2*atan(ew*ttd2*sqrt(1+(m^2)/(p^2)*(sinh(w))/(cosh(w)-cospi(t)*sinh(w))))/pi

    else # no approximation

        pB = sqrt((cosh(w)*sqrt(p^2+m^2)-sinh(w)*p*cospi(t))^2-m^2)

        if pB < 0e0
            println("pB < 0e0, $pB, setting pB=0e0")
            pB = 0e0
        end

        if w <=-4e0
            # keep to first order
            tB = 2*atan(exp(w)*p*sinpi(t),sqrt(p^2+m^2)+p*cospi(t))/pi

        elseif w >= 4e0
            # keep to first order 
            tB = 2*atan(exp(w)*(sqrt(p^2+m^2)-p*cospi(t)),p*sinpi(t))/pi 

        else # no approximation

            y = pB+exp(w)*(sqrt(p^2+m^2)-p*cospi(t))/2-exp(-w)*(p*cospi(t)+sqrt(p^2+m^2))/2
            x = pB-exp(w)*(sqrt(p^2+m^2)-p*cospi(t))/2+exp(-w)*(p*cospi(t)+sqrt(p^2+m^2))/2

            if y<0e0
                println("")
                println("warning: x<0,=$x, setting x=0e0")
                println("x=$x")
                println("y=$y")
                println("pB=$pB")
                println("w=$w")
                println("p=$p")
                println("t=$t")
                println("")
                y = 0e0
            end

            # acos(x) can lead to a domain error with x > 1 due to floating point precision, atan approach better 
            #tB = acos((cosh(w)*p*cospi(t)-sinh(w)*sqrt(p^2+m^2))/pB)/pi

            #tBtest = 2*atan(sqrt(-p*cosh(w)*cospi(t)+sqrt(p^2+m^2)*sinh(w)+pB),sqrt(+p*cosh(w)*cospi(t)-sqrt(p^2+m^2)*sinh(w)+pB))/pi

            tB = 2*atan(sqrt(y),sqrt(x))/pi

        end

    end

    pBv[1] = pB
    pBv[2] = cospi(tB)
    pBv[3] = hB
    pBv[4] = tB

    if isnan(pBv[1])
        println("pB=$pB")
        println("w=$w")
        println("p=$p")
        println("t=$t")
    end

end

"""
    DopplerFactor!(pv::Vector{Float64},m::Float64,w::Float64)

Computes the doppler factor aka cosine angle jacobian for a particle of mass `m` between the CM frame and lab frame, for a boost with rapididty `w` in the z direction. The result is stored in the fifth entry of the vector `pv`. Which is expected to have components [p,cos(theta),phi,theta,doppler].
"""
function DopplerFactor!(pv::Vector{Float64},m::Float64,w::Float64)

    p::Float64 = pv[1]
    t::Float64 = pv[4]

    if m == 0e0

        ew = exp(w)

        DF = 1/(ew*sinpi(t/2)^2+cospi(t/2)^2/ew)^2

    elseif p/m <= 1e-4 # small p approximation 
        # valid to at least second order in p/m

        DF = p^2/m^2 * cospi(t)/sinh(w)^2

    elseif p/m >= 1e4 # large p approximation
        # valid to at least second order in m/p

        ew = exp(w)
        DF = 1/(ew*sinpi(t/2)^2+cospi(t/2)^2/ew)^2
        DF += DF^2*(m^2/p^2)*sinh(w)*(2*cospi(t)*cosh(w)+(sinpi(t)-4)*sinh(w))/2

    else # no approximation

        DF = p^2*(-p*cosh(w)+cospi(t)*sqrt(p^2+m^2)*sinh(w))*(-m^2+(sqrt(p^2+m^2)*cosh(w)-p*cospi(t)*sinh(w))^2)^(-3/2)

    end

    pv[5] = abs(DF)

end

eng_error = -4.440892098500626e-16
mom_z_error = 0.00039123319021949765
mom_x_error = 1.3199434366839569e-5
w=4.0244559713931505
p1 = 128831.68341951884
ct1 = 0.23036418495849276
h1 = 0.18025009590902452
p2 = 46.55860148874508
ct2 = -0.30034740299314544
h2 = 1.4088925278594582
p3 = 41.15085565200894
ct3 = -0.5129989373463408
h3 = 1.025737771716258
p4 = 128837.07902061428
ct4 = 0.23031977234590043
h4 = 0.18017578422786484

p1*ct1
p2*ct2
p3*ct3
p4*ct4

p1*ct1+p2*ct2
p3*ct3+p4*ct4

sqrt(p1^2+1^2)+sqrt(p2^2)
sqrt(p3^2+1^2)+sqrt(p4^2)

p1*sin(acos(ct1))*cospi(h1)+p2*sin(acos(ct2))*cospi(h2)
p3*sin(acos(ct3))*cospi(h3)+p4*sin(acos(ct4))*cospi(h4)

