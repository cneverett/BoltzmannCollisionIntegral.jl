"""
    SyncEvaluateMultiThread!(userSyncInputMultiThread)

Function to run the Monte Carlo integration of the S array in a serial enviroment. 
"""
function SyncEvaluateMultiThread(userInputSyncMultiThread::Tuple{Tuple{String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64, Float64},Int64,Int64,Int64,String,String})

    # ======= Load User Parameters ======= #

    (Parameters,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName) = userInputSyncMultiThread;
    (name1,name2,mu1,mu2,z1,z2,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,BMag) = Parameters;

    # ==================================== #

    # ========= Valid Grids? =============== #
        
        if u1_grid == "b" && u1_num%2 == 0 || u2_grid == "b" && u2_num%2 == 0
            error("Binary grid must have odd number of bins")
        end

    # ====================================== #

    # ======== Load/Create Files ========= #
        
        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            SAtotal = f["STotal"];
            SAtally = f["STally"];
            SMatrix = f["SMatrix"];
            #p1Max = f["p3Max"];
            #t1MinMax = f["t3MinMax"];
            close(f)
        else
            SAtotal = zeros(Float64,p1_num,u1_num,p2_num,u2_num);  
            SAtally = zeros(UInt32,p1_num,u1_num,p2_num,u2_num);
            SMatrix = zeros(Float64,p1_num,u1_num,p2_num,u2_num);
            #pMax = zeros(Float64,p1_num,u1_num,p2_num,u2_num);
            #tMinMax = zeros(Float64,2,p1_num,u1_num,p2_num,u2_num);
            #fill!(@view(tMinMax[1,:,:,:,:]),1e0);
        end

    # ================================= #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p2_num]    

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        workers  = [SyncMonteCarloAxi_MultiThread!(SAtotal,SAtally,#=pMax,tMinMax,=#ArrayOfLocks,Parameters,numTiterPerThread,numSiterPerThread,nThreads) for _ in 1:nThreads]

        wait.(workers) # Allow all workers to finish

    # ===================================== #


    # ===== Calculate S and T Matrices === #

        # preallocate
        SMatrixOld = SMatrix;
        fill!(SMatrix,0e0);
        
        # divide element wise by tallys
        @. SMatrix = SAtotal / SAtally
        replace!(SMatrix,NaN=>0e0); # remove NaN caused by /0e0

        # Angle / Momentum Ranges
        u1val = bounds(u_low,u_up,u1_num,u1_grid)
        u2val = bounds(u_low,u_up,u2_num,u2_grid)
        p1val = bounds(p1_low,p1_up,p1_num,p1_grid)
        p2val = bounds(p2_low,p2_up,p2_num,p2_grid)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactorsSync1!(SMatrix,p1val,u1val,p2val,u2val)      # applies phase space factors for symmetries                  
        SyncSymmetry!(SMatrix)   # initial states are symmetric -> apply symmetry of interaction to improve MC values
        PhaseSpaceFactorsSync2!(SMatrix,p1val,u1val)            # corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

        # output a measure of convergence, i.e. new-old/old
        SConverge = (SMatrix .- SMatrixOld)./SMatrixOld

    # ===================================== #

    # ========== Save Arrays ============== #
            
        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal",SAtotal)
        write(f,"STally",SAtally)
        write(f,"SMatrix",SMatrix)
        #write(f,"pMax",pMax)
        #write(f,"tMinMax",tMinMax)
        write(f,"SConverge",SConverge)
        write(f,"Parameters",Parameters)
        close(f)

    # ===================================== #

end # function 