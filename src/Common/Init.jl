#= 
user defined constant parameters for STIntegration
=#

#include("ParticleData.jl")

# ----------- Particle selection ---------------- # 

    const name1::String = "Sph";
    const name2::String = "Sph";
    const name3::String = "Sph";
    const name4::String = "Sph";

# ---------------------------------------------- #

# ---------- Integration Parameters ------------ #

    # integration time approx 250ns per itteration 

    # For Serial 
    numTiter::Int64 = 1000000;    # number of T matrix itterations i.e. random p1 p2 points
    numSiter::Int64 = 100;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point
    
    #For MultiThread
    numTiterPerThread::Int64 = 10;    # number of T matrix itterations i.e. random p1 p2 points
    numSiterPerThread::Int64 = 100;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point  
    nThreads::Int64 = 1;   

# ---------------------------------------------- #

##################################################

# ----- DO NOT EDIT THIS SECTION --------------- #

    # Set Masses
    const mu1::Float32 = eval(Symbol(name1*"Data")).mu
    const mu2::Float32 = eval(Symbol(name2*"Data")).mu
    const mu3::Float32 = eval(Symbol(name3*"Data")).mu
    const mu4::Float32 = eval(Symbol(name4*"Data")).mu


    # Set Momentum Space Grids
    
    #const log10pspace::Bool = true;

    #if (log10pspace == true)
        # momentum space grids are defined in log10 space
        const p3l::Float32 = eval(Symbol(name3*"Data")).pl
        const p3u::Float32 = eval(Symbol(name3*"Data")).pu
        const nump3::Int64 = eval(Symbol(name3*"Data")).nump

        const p1l::Float32 = eval(Symbol(name1*"Data")).pl
        const p1u::Float32 = eval(Symbol(name1*"Data")).pu
        const nump1::Int64 = eval(Symbol(name1*"Data")).nump

        const p2l::Float32 = eval(Symbol(name2*"Data")).pl
        const p2u::Float32 = eval(Symbol(name2*"Data")).pu
        const nump2::Int64 = eval(Symbol(name2*"Data")).nump

    #elseif (log10pspace == false)
        # momentum space grids are defined in normal space
    #    const p3u::Float32 = 1f4;
    #    const p3l::Float32 = 1f-2;
    #    const nump3::UInt32 = UInt32(6);

    #    const p1u::Float32 = 1f4;
    #    const p1l::Flaot32 = 1f-2;
    #    const nump1::UInt32 = UInt32(6);

    #    const p2u::Flaot32 = 1f4;
    #    const p2l::Float32 = 1f-2;
    #    const nump2::UInt32 = UInt32(6);
    #else 
    #    error("log10pspace not defined")
    #end

    # cos(theta) 
    const t3l::Float32 = eval(Symbol(name3*"Data")).tl
    const t3u::Float32 = eval(Symbol(name3*"Data")).tu
    const numt3::Int64 = eval(Symbol(name3*"Data")).numt

    const t1l::Float32 = eval(Symbol(name1*"Data")).tl
    const t1u::Float32 = eval(Symbol(name1*"Data")).tu
    const numt1::Int64 = eval(Symbol(name1*"Data")).numt
    
    const t2l::Float32 = eval(Symbol(name2*"Data")).tl
    const t2u::Float32 = eval(Symbol(name2*"Data")).tu
    const numt2::Int64 = eval(Symbol(name2*"Data")).numt

# ---------------------------------------------- #

##################################################

# --------------- File Location ---------------- #

    # file name contains all the information of discretisation in form pl3#pu3#nump3#p1l#p1u#nump1#p2l#p2u#nump2#numt3#numt1#numt2 (tu/tl not needed as bound is always [-1,1])

    fileLocation = pwd()*"\\Data"
    fileName = name1*name2*name3*name4*"#"*string(Int(p3l))*"#"*string(Int(p3u))*"#"*string(nump3)*"#"*string(Int(p1l))*"#"*string(Int(p1u))*"#"*string(nump1)*"#"*string(Int(p2l))*"#"*string(Int(p2u))*"#"*string(nump2)*"#"*string(numt3)*"#"*string(numt1)*"#"*string(numt2)*".jld2"

# ---------------------------------------------- #