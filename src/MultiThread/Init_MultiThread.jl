# ----------- Particle selection ---------------- # 

    const name1::String = "Sph";
    const name2::String = "Sph";
    const name3::String = "Sph";
    const name4::String = "Sph";

# ---------------------------------------------- #

# ---------- Integration Parameters ------------ #

# integration time approx 250ns per itteration 

    const numTiterPerThread::Int64 = 10;    # number of T matrix itterations i.e. random p1 p2 points
    const numSiterPerThread::Int64 = 100;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point  
    const nThreads::Int64 = 1;   

# ---------------------------------------------- #

# --------------- File Location ---------------- #

    fileLocation = pwd()*"\\Data"
    fileName = filename = name1*name2*name3*name4*".jld2"

# ---------------------------------------------- #