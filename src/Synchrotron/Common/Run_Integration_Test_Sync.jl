#=
This script is provided as an example of how to run the BinaryInteractionSpectra.jl Synchrotron module.
=#

# First select the synchrotron emitting particle (particle 2) you want to evaluate the synchrotron emissions. 

    name2::String = "Ele";

# Define the Momentum space discretisations for particle 2 and particle 1 (particle 1 is the emitted photon). This includes the upper and lower momentum bounds (log10) and the number of bins for momentum magnitude and cos(theta) and must be of the format `pl_name`, pl_name`, `nump_name` and `numt_name` where `name` is the abreviated three letter name of the particle species. 

    pl_Ele = -5e0
    pu_Ele = 4e0
    nump_Ele = Int64(8*(pu_Ele-pl_Ele))
    numt_Ele = 8

    pl_Pho = -14e0
    pu_Pho = 4e0
    nump_Pho = Int64(4*(pu_Pho-pl_Pho))
    numt_Pho = 8

# Define the Magnetic Field Strength (in Tesla)  

    BMag::Float64 = 1e-6;

# ==== DO NOT EDIT THIS SECTION ======= #
# ===================================== #

    # This section takes all the information above, and geneartes a tuple of the parameters to be passed to the evaluation functions.

    name1::String = "Pho";
    p1l::Float64 = getfield(Main,Symbol("pl_"*name1))
    p1u::Float64 = getfield(Main,Symbol("pu_"*name1))
    nump1::Int64 = getfield(Main,Symbol("nump_"*name1))
    numt1::Int64 = getfield(Main,Symbol("numt_"*name1))
    p2l::Float64 = getfield(Main,Symbol("pl_"*name2))
    p2u::Float64 = getfield(Main,Symbol("pu_"*name2))
    nump2::Int64 = getfield(Main,Symbol("nump_"*name2))
    numt2::Int64 = getfield(Main,Symbol("numt_"*name2))

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
# By default the file name contains all the information of discretisation in form p1l#p1u#nump1#p2l#p2u#nump2#numt1#numt2, and it to be stored in a folder called "Data" in the current working directory.

    fileLocation::String = pwd()*"\\Data";
    fileName::String = "sync"*name2*"#"*string(p1l)*"#"*string(p1u)*"#"*string(nump1)*"#"*string(p2l)*"#"*string(p2u)*"#"*string(nump2)*"#"*string(numt1)*"#"*string(numt2)*".jld2";


# Now run the evaluation functions. (Comment out as needed)

    # ===== DO NOT EDIT THIS SECTION ===== #
    
        userInputSyncSerial = (name2,p1l,p1u,nump1,numt1,p2l,p2u,nump2,numt2,numTiter,numSiter,fileLocation,fileName,BMag);
        userInputSyncMultiThread = (name2,p1l,p1u,nump1,numt1,p2l,p2u,nump2,numt2,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName,BMag);

    # ===================================== #

    # For Serial Evaluation
    #@time SyncEvaluateSerial(userInputSyncSerial)

    # For MultiThread Evaluation
    @time SyncEvaluateMultiThread(userInputSyncMultiThread)

# Once evaluation is complete all the data stored in the output file can be loaded into the user workspace using the following function:

    #(Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName);


