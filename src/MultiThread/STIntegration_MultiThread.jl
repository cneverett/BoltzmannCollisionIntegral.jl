#= Script for running the ST integration and returning data arrays =#

#=
    Optimisation of the multi-threaded version would not have been possible without the kind guidence of those on the Julia Forums: https://discourse.julialang.org/t/fast-multi-threaded-array-editing-without-data-races/114863/41
    In particular the assistance of users: mbauman, adienes, Oscar_Smith, Satvik, Salmon, sgaure and foobar_lv2
=#

"""
    SpectraEvaluateMultiThread(userInputMultiThread)

Function to run the Monte Carlo integration of the S and T arrays in a multi-threaded environment. The function will run the Monte Carlo integration in parallel across the number of threads specified in the global variable nThreads. The function will then calculate the S and T matrices and save the results to a file.
"""
function SpectraEvaluateMultiThread(userInputMultiThread::Tuple{String,String,String,String,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,String,String})

    # ========= Load user Parameters ======= #

    (name1,name2,name3,name4,p3_low,p3_up,num_p3,p4_low,p4_up,num_p4,p1_low,p1_up,num_p1,p2_low,p2_up,num_p2,num_u3,num_u4,num_u1,num_u2,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName) = userInputMultiThread

    # ====================================== #

    # ===== Are states Distinguishable ===== #

        Indistinguishable_12::Bool = name1 == name2
        Indistinguishable_34::Bool = name3 == name4

    # ====================================== #

    # ========= Load/Create Files ========== #

        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            SAtotal3 = f["STotal3"];
            SAtally3 = f["STally3"];
            SMatrix3 = f["SMatrix3"];
            p3Max = f["p3Max"];
            u3MinMax = f["u3MinMax"];
            SAtotal4 = f["STotal4"];
            SAtally4 = f["STally4"];
            SMatrix4 = f["SMatrix4"];
            p4Max = f["p4Max"];
            u4MinMax = f["u4MinMax"];
            TAtotal = f["TTotal"];
            TAtally = f["TTally"];
            TMatrix1 = f["TMatrix1"];
            TMatrix2 = f["TMatrix2"];
            close(f)
        else
            SAtotal3 = zeros(Float64,(num_p3+1),num_u3,num_p1,num_u1,num_p2,num_u2); 
            SAtotal4 = zeros(Float64,(num_p4+1),num_u4,num_p1,num_u1,num_p2,num_u2); 
            TAtotal = zeros(Float64,num_p1,num_u1,num_p2,num_u2);
            SAtally3 = zeros(UInt32,num_u3,num_p1,num_u1,num_p2,num_u2);
            SAtally4 = zeros(UInt32,num_u4,num_p1,num_u1,num_p2,num_u2)
            TAtally = zeros(UInt32,num_p1,num_u1,num_p2,num_u2);
            SMatrix3 = zeros(Float64,(num_p3+1),num_u3,num_p1,num_u1,num_p2,num_u2);
            SMatrix4 = zeros(Float64,(num_p4+1),num_u4,num_p1,num_u1,num_p2,num_u2);
            TMatrix1 = zeros(Float64,num_p1,num_u1,num_p2,num_u2);
            TMatrix2 = zeros(Float64,num_p2,num_u2,num_p1,num_u1);
            p3Max = zeros(Float64,num_u3,num_p1,num_u1,num_p2,num_u2);
            u3MinMax = zeros(Float64,2,(num_p3+1),num_p1,num_u1,num_p2,num_u2);
            p4Max = zeros(Float64,num_u4,num_p1,num_u1,num_p2,num_u2);
            u4MinMax = zeros(Float64,2,(num_p4+1),num_p1,num_u1,num_p2,num_u2);
            fill!(@view(u3MinMax[1,:,:,:,:,:]),1e0);
            fill!(@view(u3MinMax[2,:,:,:,:,:]),-1e0);
            fill!(@view(u4MinMax[1,:,:,:,:,:]),1e0);
            fill!(@view(u4MinMax[2,:,:,:,:,:]),-1e0);
        end

    # ====================================== #

    # ======= Define Cross Section Functions Based on Particle Selections ========= #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma = getfield(BoltzmannCollisionIntegral,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt = getfield(BoltzmannCollisionIntegral,name_dsigmadt)

    # ============================================================================ #

    # ===== Set Particle (normalised) Masses) and Parameters ====== #

        mu1::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
        mu2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
        mu3::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
        mu4::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

        Parameters = (mu1,mu2,mu3,mu4,p3_low,p3_up,num_p3,p4_low,p4_up,num_p4,p1_low,p1_up,num_p1,p2_low,p2_up,num_p2,num_u3,num_u4,num_u1,num_u2)

    # ============================================================= #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:num_p1]    

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        # Set up workers
        workers = [STMonteCarloAxi_MultiThread!(SAtotal3,SAtotal4,TAtotal,SAtally3,SAtally4,TAtally,ArrayOfLocks,p3Max,p4Max,u3MinMax,u4MinMax,sigma,dsigmadt,Parameters,numTiterPerThread,numSiterPerThread) for _ in 1:nThreads]
        
        wait.(workers) # Allow all workers to finish
   
    # ===================================== #

    # ===== Calculate S and T Matricies === #

        #= these operations are run in serial =#

        # preallocate
        SMatrixOldSum3 = dropdims(sum(SMatrix3,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix3,0e0);
        SMatrixOldSum4 = dropdims(sum(SMatrix4,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix4,0e0);
        TMatrixOldSum = dropdims(sum(TMatrix1,dims=(3,4)),dims=(3,4));
        fill!(TMatrix1,0e0);

        # divide element wise by tallys
        if Indistinguishable_34 == true # 34 are identical
            @. SAtally3 = SAtally3 + SAtally4
            @. SAtotal3 = SAtotal3 + SAtotal4
            @. p3Max = max(p3Max,p4Max)
            @view(u3MinMax[1,:,:,:,:,:]) .= min.(u3MinMax[1,:,:,:,:,:],u4MinMax[1,:,:,:,:,:])
            @view(u3MinMax[2,:,:,:,:,:]) .= max.(u3MinMax[2,:,:,:,:,:],u4MinMax[2,:,:,:,:,:])
            # reset arrays to avoid overcounting when multiple runs are made
            fill!(SAtally4,UInt32(0))
            fill!(SAtotal4,0e0)
            fill!(p4Max,0e0)
            fill!(@view(u4MinMax[1,:,:,:,:,:]),1e0)
            fill!(@view(u4MinMax[2,:,:,:,:,:]),-1e0)
        end
        
        if (Indistinguishable_34 == false) && (mu3 == mu4)  # system symmetric in 34 interchange but particles not identical
            @. SAtally3 = SAtally3 + SAtally4
            @. SAtally4 = SAtally3
            @. SAtotal3 = SAtotal3 + SAtotal4
            @. SAtotal4 = SAtotal3
            @. p3Max = max(p3Max,p4Max)
            @. p4Max = p3Max
            @view(u3MinMax[1,:,:,:,:,:]) .= min.(u3MinMax[1,:,:,:,:,:],u4MinMax[1,:,:,:,:,:])
            @view(u3MinMax[2,:,:,:,:,:]) .= max.(u3MinMax[2,:,:,:,:,:],u4MinMax[2,:,:,:,:,:])
            @. u4MinMax = u3MinMax
        end

        for i in axes(SMatrix3,1)
            @. @view(SMatrix3[i,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:]) / SAtally3
        end
        replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
        for i in axes(SMatrix4,1)
            @. @view(SMatrix4[i,:,:,:,:,:]) = @view(SAtotal4[i,:,:,:,:,:]) / SAtally4
        end
        replace!(SMatrix4,NaN=>0e0); # remove NaN caused by /0e0
    
        TMatrix1 = TAtotal ./ TAtally;
        replace!(TMatrix1,NaN=>0e0);

        # Angle / Momentum Grid Bounds
        u3_bounds = bounds(u_up,u_low,num_u3,"u")
        u4_bounds = bounds(u_up,u_low,num_u4,"u")
        u1_bounds = bounds(u_up,u_low,num_u1,"u")
        u2_bounds = bounds(u_up,u_low,num_u2,"u")
        p3_bounds = bounds(p3_up,p3_low,num_p3,"l")
        p4_bounds = bounds(p4_up,p4_low,num_p4,"l")
        p1_bounds = bounds(p1_up,p1_low,num_p1,"l")
        p2_bounds = bounds(p2_up,p2_low,num_p2,"l")

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix3,SMatrix4,TMatrix1,u3_bounds,u4_bounds,p1_bounds,u1_bounds,p2_bounds,u2_bounds,Indistinguishable_12)      # applies phase space factors for symmetries                  
        STSymmetry!(SMatrix3,SMatrix4,TMatrix1,mu1,mu2)   # initial states are symmetric -> apply symmetry of interaction to improve MC values
        if Indistinguishable_12 == false
            perm = [3,4,1,2]
            TMatrix2 = permutedims(TMatrix1,perm)
        end
        PhaseSpaceFactors2!(SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3_bounds,u3_bounds,p4_bounds,u4_bounds,p1_bounds,u1_bounds,p2_bounds,u2_bounds)            # corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

        # output a measure of convergence, i.e. new-old/old
        SMatrixSum3 = dropdims(sum(SMatrix3,dims=(3,4,5,6)),dims=(3,4,5,6));
        SConverge3 = (SMatrixSum3 .- SMatrixOldSum3)./SMatrixOldSum3
        SMatrixSum4 = dropdims(sum(SMatrix4,dims=(3,4,5,6)),dims=(3,4,5,6));
        SConverge4 = (SMatrixSum4 .- SMatrixOldSum4)./SMatrixOldSum4
        TMatrixSum = dropdims(sum(TMatrix1,dims=(3,4)),dims=(3,4));
        TConverge = (TMatrixSum .- TMatrixOldSum)./TMatrixOldSum

    # ===================================== # 

    # ========== Save Arrays ============== #

        OutputParameters = (name1,name2,name3,name4,p3_low,p3_up,num_p3,p4_low,p4_up,num_p4,p1_low,p1_up,num_p1,p2_low,p2_up,num_p2,num_u3,num_u4,num_u1,num_u2)
        
        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal3",SAtotal3)
        write(f,"STally3",SAtally3)
        write(f,"SMatrix3",SMatrix3)
        write(f,"p3Max",p3Max)
        write(f,"u3MinMax",u3MinMax)
        write(f,"SConverge3",SConverge3)

        write(f,"STotal4",SAtotal4)
        write(f,"STally4",SAtally4)
        write(f,"SMatrix4",SMatrix4)
        write(f,"p4Max",p4Max)
        write(f,"u4MinMax",u4MinMax)
        write(f,"SConverge4",SConverge4)

        write(f,"TTotal",TAtotal)
        write(f,"TTally",TAtally)
        write(f,"TMatrix1",TMatrix1)
        write(f,"TMatrix2",TMatrix2)

        write(f,"TConverge",TConverge)

        write(f,"Parameters",OutputParameters)
        close(f)

    # ===================================== #

    return nothing

end # function