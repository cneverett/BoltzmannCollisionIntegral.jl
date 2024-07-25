#=
This script is provided as an example of how to run the BinaryInteractionSpectra.jl module.
=#

# First select the binary interaction (12->34) you want to evaluate the collision integral for by setting the names of particles 1 through 4. The pairs of names (12) and (34) should be in alphabetical order.

    name1::String = "Sph";
    name2::String = "Sph";
    name3::String = "Sph";
    name4::String = "Sph";

# Define the Momentum space discretisations for each particle species named. This includes the upper and lower momentum bounds and the number of bins for momentum magnitude and cos(theta) and must be of the format `pl_name`, pl_name`, `nump_name` and `numt_name` where `name` is the abreviated three letter name of the particle species. 

    pl_Sph = -5f0
    pu_Sph = 4f0
    nump_Sph = 72
    numt_Sph = 8

# ==== DO NOT EDIT THIS SECTION ======= #
# ===================================== #

    # This section takes all the information above, and geneartes a tuple of the parameters to be passed to the evaluation functions.

    p1l::Float64 = getfield(Main,Symbol("pl_"*name1))
    p1u::Float64 = getfield(Main,Symbol("pu_"*name1))
    nump1::Int64 = getfield(Main,Symbol("nump_"*name1))
    numt1::Int64 = getfield(Main,Symbol("numt_"*name1))
    p2l::Float64 = getfield(Main,Symbol("pl_"*name2))
    p2u::Float64 = getfield(Main,Symbol("pu_"*name2))
    nump2::Int64 = getfield(Main,Symbol("nump_"*name2))
    numt2::Int64 = getfield(Main,Symbol("numt_"*name2))
    p3l::Float64 = getfield(Main,Symbol("pl_"*name3))
    p3u::Float64 = getfield(Main,Symbol("pu_"*name3))
    nump3::Int64 = getfield(Main,Symbol("nump_"*name3))
    numt3::Int64 = getfield(Main,Symbol("numt_"*name3))
    p4l::Float64 = getfield(Main,Symbol("pl_"*name4))
    p4u::Float64 = getfield(Main,Symbol("pu_"*name4))
    nump4::Int64 = getfield(Main,Symbol("nump_"*name4))
    numt4::Int64 = getfield(Main,Symbol("numt_"*name4))

# ===================================== #
# ===================================== #

# Define the number of Monte-Carlo samples to take. numTiter gives the number of random sets of (p1,p2) vectors that are sampled, while numSiter given the number of random (p3) vectors that are sampled for each set of (p1,p2) vectors. For good convergence, these values should be much larger than the number of bins in the momentum space discretisation. 
# As a rule of thumb, each single itteration takes approx 250ns to complete on a single core of a modern CPU

    # For Serial Evaluation
    numTiter::Int64 = nump1*numt1*nump2*numt2*1;    # number of T matrix itterations i.e. random p1 p2 points
    numSiter::Int64 = nump3*numt3*1;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point

    # For MultiThread Evaluation
    numTiterPerThread::Int64 = nump1*numt1*nump2*numt2*1;    # number of T matrix itterations i.e. random p1 p2 points. Should be > nump1*numt1*nump2*numt2 to ensure good sampling.
    numSiterPerThread::Int64 = nump3*numt3*1;        # number of S matrix iteration per T matrix iteration i.e. random p3 directions per p1 p2 point. Should be > nump3*numt3 to ensure good sampling.

# Define the number of Threads (more specifically, workers) used for MultiThreaded evaluation. (see https://docs.julialang.org/en/v1/manual/multi-threading/ for how to set up multi-threading in Julia)

    nThreads::Int64 = 10;

# Define the location where the output file (.jld2 format) is to be stored and its name.
# By default the file name contains all the information of discretisation in form pl3#pu3#nump3#p1l#p1u#nump1#p2l#p2u#nump2#numt3#numt1#numt2, and it to be stored in a folder called "Data" in the current working directory.

    fileLocation::String = pwd()*"\\Data";
    fileName::String = name1*name2*name3*name4*"#"*string(Int(p3l))*"#"*string(Int(p3u))*"#"*string(nump3)*"#"*string(Int(p1l))*"#"*string(Int(p1u))*"#"*string(nump1)*"#"*string(Int(p2l))*"#"*string(Int(p2u))*"#"*string(nump2)*"#"*string(numt3)*"#"*string(numt1)*"#"*string(numt2)*".jld2"

# Now run the evaluation functions.

    # ===== DO NOT EDIT THIS SECTION ===== #
    
        userInputSerial = (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2,numTiter,numSiter,fileLocation,fileName)
        userInputMultiThread = (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName)

    # ===================================== #

    # For Serial Evaluation
    SpectraEvaluateSerial(userInputSerial)

    # For MultiThread Evaluation
    SpectraEvaluateMultiThread(userInputMultiThread)

# Once evaluation is complete all the data stored in the output file can be loaded into the user workspace using the following function:

    (Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName);


