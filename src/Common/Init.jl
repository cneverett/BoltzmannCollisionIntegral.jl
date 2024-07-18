#= 
user defined constant parameters for STIntegration
=#

#include("ParticleData.jl")

# ----------- Particle selection ---------------- # 

    name1::String = "Sph";
    name2::String = "Sph";
    name3::String = "Sph";
    name4::String = "Sph";

# ---------------------------------------------- #

##################################################

# ----- DO NOT EDIT THIS SECTION --------------- #

    # Set (Normalised) Masses
    const mu1::Float32 = getfield(BinaryInteractionSpectra,Symbol("mu"*name1))
    const mu2::Float32 = getfield(BinaryInteractionSpectra,Symbol("mu"*name2))
    const mu3::Float32 = getfield(BinaryInteractionSpectra,Symbol("mu"*name3))
    const mu4::Float32 = getfield(BinaryInteractionSpectra,Symbol("mu"*name4))

    # Set Momentum Space Grids defined in log10 space
    const p3l::Float32 = -5f0
    const p3u::Float32 = 4f0
    const nump3::Int64 = 72

    const p1l::Float32 = -5f0
    const p1u::Float32 = 4f0
    const nump1::Int64 = 72

    const p2l::Float32 = -5f0
    const p2u::Float32 = 4f0
    const nump2::Int64 = 72

    const numt3::Int64 = 8
    const numt1::Int64 = 8
    const numt2::Int64 = 8

# ---------------------------------------------- #

##################################################

# ---------- Integration Parameters ------------ #

    # integration time approx 250ns per itteration 

    # For Serial 
    numTiter::Int64 = nump1*numt1*nump2*numt2*1;    # number of T matrix itterations i.e. random p1 p2 points
    numSiter::Int64 = nump3*numt3*1;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point

    #For MultiThread
    numTiterPerThread::Int64 = nump1*numt1*nump2*numt2*35;    # number of T matrix itterations i.e. random p1 p2 points. Should be > nump1*numt1*nump2*numt2 to ensure good sampling.
    numSiterPerThread::Int64 = nump3*numt3*35;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point. Should be > nump3*numt3 to ensure good sampling.
    nThreads::Int64 = 10;    # number of threads
    
# ---------------------------------------------- #

# --------------- File Location ---------------- #

    # file name contains all the information of discretisation in form pl3#pu3#nump3#p1l#p1u#nump1#p2l#p2u#nump2#numt3#numt1#numt2 (tu/tl not needed as bound is always [-1,1])

    fileLocation = pwd()*"\\Data"
    fileName = name1*name2*name3*name4*"#"*string(Int(p3l))*"#"*string(Int(p3u))*"#"*string(nump3)*"#"*string(Int(p1l))*"#"*string(Int(p1u))*"#"*string(nump1)*"#"*string(Int(p2l))*"#"*string(Int(p2u))*"#"*string(nump2)*"#"*string(numt3)*"#"*string(numt1)*"#"*string(numt2)*"new1.jld2"

# ---------------------------------------------- #