#= Script for running the ST integration and returning data arrays =#

#=
    Optimisation of the multi-threaded version would not have been possible without the kind guidence of those on the Julia Forums: https://discourse.julialang.org/t/fast-multi-threaded-array-editing-without-data-races/114863/41
    In particular the assistance of users: mbauman, adienes, Oscar_Smith, Satvik, Salmon, sgaure and foobar_lv2
=#

"""
    SpectraEvaluateMultiThread(userInputMultiThread)

Function to run the Monte Carlo integration of the S and T arrays in a multi-threaded environment. The function will run the Monte Carlo integration in parallel across the number of threads specified in the global variable nThreads. The function will then calculate the S and T matrices and save the results to a file.
"""
function SpectraEvaluateMultiThread(userInputMultiThread::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64},Int64,Int64,Int64,String,String})

    # ========= Load user Parameters ======= #

    (Parameters,numTiterPerThread,numSiterPerThread,nThreads,fileLocation,fileName) = userInputMultiThread
    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num) = Parameters

    # ====================================== #

    # ========= Valid Grids? =============== #
        
        if u1_grid == "b" && u1_num%2 == 0 || u2_grid == "b" && u2_num%2 == 0 || u3_grid == "b" && u3_num%2 == 0 || u4_grid == "b" && u4_num%2 == 0
            error("Binary grid must have odd number of bins")
        end

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
            SAtotal3 = zeros(Float64,(p3_num+1),u3_num,p1_num,u1_num,p2_num,u2_num); 
            SAtotal4 = zeros(Float64,(p4_num+1),u4_num,p1_num,u1_num,p2_num,u2_num); 
            TAtotal = zeros(Float64,p1_num,u1_num,p2_num,u2_num);
            SAtally3 = zeros(UInt32,u3_num,p1_num,u1_num,p2_num,u2_num);
            SAtally4 = zeros(UInt32,u4_num,p1_num,u1_num,p2_num,u2_num)
            TAtally = zeros(UInt32,p1_num,u1_num,p2_num,u2_num);
            SMatrix3 = zeros(Float64,(p3_num+1),u3_num,p1_num,u1_num,p2_num,u2_num);
            TMatrix1 = zeros(Float64,p1_num,u1_num,p2_num,u2_num);
            TMatrix2 = zeros(Float64,p2_num,u2_num,p1_num,u1_num);
            p3Max = zeros(Float64,u3_num,p1_num,u1_num,p2_num,u2_num);
            u3MinMax = zeros(Float64,2,(p3_num+1),p1_num,u1_num,p2_num,u2_num);
            SMatrix4 = zeros(Float64,(p4_num+1),u4_num,p1_num,u1_num,p2_num,u2_num);
            p4Max = zeros(Float64,u4_num,p1_num,u1_num,p2_num,u2_num);
            u4MinMax = zeros(Float64,2,(p4_num+1),p1_num,u1_num,p2_num,u2_num);
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

        #=mu1::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
        mu2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
        mu3::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
        mu4::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

        Parameters = (mu1,mu2,mu3,mu4,p3_low,p3_up,num_p3,p4_low,p4_up,num_p4,p1_low,p1_up,p1_num,p2_low,p2_up,p2_num,num_u3,num_u4,u1_num,num_u2) =#

    # ============================================================= #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:p1_num]    

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        # Set up workers
        workers = [STMonteCarloAxi_MultiThread!(SAtotal3,SAtotal4,TAtotal,SAtally3,SAtally4,TAtally,ArrayOfLocks,p3Max,p4Max,u3MinMax,u4MinMax,sigma,dsigmadt,Parameters,numTiterPerThread,numSiterPerThread) for _ in 1:nThreads]
        
        wait.(workers) # Allow all workers to finish
   
    # ===================================== #

    # ===== Calculate S and T Matrices === #

        # preallocate
        SMatrixOldSum3 = dropdims(sum(SMatrix3,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix3,0e0);
        SMatrixOldSum4 = dropdims(sum(SMatrix4,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix4,0e0);
        TMatrixOldSum = dropdims(sum(TMatrix1,dims=(3,4)),dims=(3,4));
        fill!(TMatrix1,0e0);

        # divide element wise by tallys
        if Indistinguishable_34 == true
            @. SAtally3 = SAtally3 + SAtally4
            @. SAtotal3 = SAtotal3 + SAtotal4
            for i in axes(SMatrix3,1)
                @. @view(SMatrix3[i,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:]) / SAtally3
            end
            replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            @. p3Max = max(p3Max,p4Max)
            @view(u3MinMax[1,:,:,:,:,:]) .= min.(u3MinMax[1,:,:,:,:,:],u4MinMax[1,:,:,:,:,:])
            @view(u3MinMax[2,:,:,:,:,:]) .= max.(u3MinMax[2,:,:,:,:,:],u4MinMax[2,:,:,:,:,:])
            # reset arrays to avoid overcounting when multiple runs are made
                fill!(SAtally4,UInt32(0))
                fill!(SAtotal4,0e0)
                fill!(p4Max,0e0)
                fill!(@view(u4MinMax[1,:,:,:,:,:]),1e0);
                fill!(@view(u4MinMax[2,:,:,:,:,:]),-1e0);
        elseif mu3 == mu4 # system symmetric in 34 interchange
            @. SAtally3 = SAtally3 + SAtally4
            @. SAtally4 = SAtally3
            @. SAtotal3 = SAtotal3 + SAtotal4
            @. SAtotal4 = SAtotal3
            for i in axes(SMatrix3,1)
                @. @view(SMatrix3[i,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:]) / SAtally3
            end
            replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            @. SMatrix4 = SMatrix3
            @. p3Max = max(p3Max,p4Max)
            @. p4Max = p3Max
            @view(u3MinMax[1,:,:,:,:,:]) .= min.(u3MinMax[1,:,:,:,:,:],u4MinMax[1,:,:,:,:,:])
            @view(u3MinMax[2,:,:,:,:,:]) .= max.(u3MinMax[2,:,:,:,:,:],u4MinMax[2,:,:,:,:,:])
            @. u4MinMax = u3MinMax
        else
            for i in axes(SMatrix3,1)
                @. @view(SMatrix3[i,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:]) / SAtally3
            end
            replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            for i in axes(SMatrix4,1)
                @. @view(SMatrix4[i,:,:,:,:,:]) = @view(SAtotal4[i,:,:,:,:,:]) / SAtally4
            end
            replace!(SMatrix4,NaN=>0e0); # remove NaN caused by /0e0
        end
        TMatrix1 = TAtotal ./ TAtally;
        replace!(TMatrix1,NaN=>0e0);

        # Angle / Momentum Ranges
        u3val = bounds(u_up,u_low,u3_num,u3_grid)
        u4val = bounds(u_up,u_low,u4_num,u4_grid)
        u1val = bounds(u_up,u_low,u1_num,u1_grid)
        u2val = bounds(u_up,u_low,u2_num,u2_grid)
        p3val = bounds(p3_low,p3_up,p3_num,p3_grid)
        p4val = bounds(p4_low,p4_up,p4_num,p4_grid)
        p1val = bounds(p1_low,p1_up,p1_num,p1_grid)
        p2val = bounds(p2_low,p2_up,p2_num,p2_grid)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix3,SMatrix4,TMatrix1,u3val,u4val,p1val,u1val,p2val,u2val,Indistinguishable_12)      # applies phase space factors for symmetries                  
        STSymmetry!(SMatrix3,SMatrix4,TMatrix1,mu1,mu2)   # initial states are symmetric -> apply symmetry of interaction to improve MC values
        if Indistinguishable_12 == false
            perm = [3,4,1,2]
            TMatrix2 = permutedims(TMatrix1,perm)
        end
        PhaseSpaceFactors2!(SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3val,u3val,p4val,u4val,p1val,u1val,p2val,u2val)            # corrects phase space factors for application in kinetic models
                                            
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

        write(f,"Parameters",Parameters)
        close(f)

    # ===================================== #


    return nothing

end # function