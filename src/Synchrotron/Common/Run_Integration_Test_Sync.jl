#=
This script is provided as an example of how to run the BinaryInteractionSpectra.jl Synchrotron module.
=#

# First select the synchrotron emitting particle (particle 2) you want to evaluate the synchrotron emissions. 

    name2::String = "Ele";

# Define the Momentum space discretisation for particle 2 and particle 1 (particle 1 is the emitted photon).  This includes the upper and lower momentum bounds, grid type "l", "u", or "b", and the number of bins. Must be of the format `p_low_name`, p_low_name`, `p_grid_name`, `p_num_name`, `u_grid_name` and `u_num_name` where `name` is the abbreviated three letter name of the particle species. 

    p_low_Ele = -5e0
    p_up_Ele = 4e0
    p_grid_Ele = "l"
    p_num_Ele = 72
    u_grid_Ele = "u"
    u_num_Ele = 8

    pl_Pho = -14e0
    pu_Pho = 4e0
    nump_Pho = 72
    numt_Pho = 8

# Define the Magnetic Field Strength (in Tesla)  

    BMag::Float64 = 1e-6;

# ==== DO NOT EDIT THIS SECTION ======= #
# ===================================== #

    # This section takes all the information above, and generates a tuple of the parameters to be passed to the evaluation functions.

    name1::String = "Pho";
    p1_low::Float64 = getfield(Main,Symbol("p_low_"*name1))
    p1_up::Float64 = getfield(Main,Symbol("p_up_"*name1))
    p1_grid::String = getfield(Main,Symbol("p_grid_"*name1))
    p1_num::Int64 = getfield(Main,Symbol("p_num_"*name1))
    u1_grid::String = getfield(Main,Symbol("u_grid_"*name1))
    u1_num::Int64 = getfield(Main,Symbol("u_num_"*name1))
    p2_low::Float64 = getfield(Main,Symbol("p_low_"*name2))
    p2_up::Float64 = getfield(Main,Symbol("p_up_"*name2))
    p2_grid::String = getfield(Main,Symbol("p_grid_"*name2))
    p2_num::Int64 = getfield(Main,Symbol("p_num_"*name2))
    u2_grid::String = getfield(Main,Symbol("u_grid_"*name2))
    u2_num::Int64 = getfield(Main,Symbol("u_num_"*name2))

    mu1::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
    mu2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))

    z1::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("z"*name1))
    z2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("z"*name2))

    Parameters = (name1,name2,mu1,mu2,z1,z2,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,BMag)

# ===================================== #
# ===================================== #

# Define the number of Monte-Carlo samples to take. numSiter gives the number of random sets of (p1,p2) vectors that are sampled. For good convergence, these values should be much larger than the number of bins in the momentum space discretisation. 

    # For Serial Evaluation
    numTiter::Int64 = nump2*numt2*10;    # number of emitting particle states to sample
    numSiter::Int64 = nump1*numt1*100;    # number of photon states per emitting particle state to sample

    # For MultiThread Evaluation
    numTiterPerThread::Int64 = nump2*numt2*100;    # number of emitting particle states to sample per thread
    numSiterPerThread::Int64 = nump1*numt1*10000;    # number of photon states per emitting particle state to sample per thread

# Define the number of Threads (more specifically, workers) used for MultiThreaded evaluation. (see https://docs.julialang.org/en/v1/manual/multi-threading/ for how to set up multi-threading in Julia)

    nThreads::Int64 = 10;

# Define the location where the output file (.jld2 format) is to be stored and its name.
# By default the file name contains all the information of discretisation in form "sync p1_low-p1_up p1_grid p1_num#...#u2_grid u2_num.jld2", and it to be stored in a folder called "Data" in the current working directory.

    fileLocation::String = pwd()*"\\Data";
    fileName::String = "sync"*name2*"#"*string(p1_low)*"-"*string(p1_up)*p1_grid*string(p1_num)*"#"*string(p2_low)*"-"*string(p2_up)*p2_grid*string(p2_num)*"#"*u1_grid*string(u1_num)*"#"*u2_grid*string(u2_num)*".jld2";

# Now run the evaluation functions. (Comment out as needed)

    # ===== DO NOT EDIT THIS SECTION ===== #
    
        userInputSyncSerial = (Parameters,numTiter,numSiter,fileLocation,fileName);
        userInputSyncMultiThread = (Parameters,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName);

    # ===================================== #

    # For Serial Evaluation
    #@time SyncEvaluateSerial(userInputSyncSerial)

    # For MultiThread Evaluation
    @time SyncEvaluateMultiThread(userInputSyncMultiThread)

# Once evaluation is complete all the data stored in the output file can be loaded into the user workspace using the following function:

    #(Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName);


