#= Script for running the ST integration and returning data arrays =#

    #include("../Common/MyPhysicalConstants.jl")
    #    include("../Common/ParticleData.jl")
    #    include("../Common\\Init.jl")
    #    include("../Common\\DifferentialCrossSectionFunctions.jl")
    #    include("../Common\\Momentum3Values.jl")
    #    include("../Common\\RandomPointMomentum.jl")
    #    include("../Common\\RandomPointSphere.jl")
    #    include("../Common/MandelstramChecks.jl")
    #    include("../Common\\STValue.jl")
    #    include("../Common/UsefulGridValueFunctions.jl")
    #    include("../Common/PhaseSpaceFactors.jl")
    #    include("../Common/Location.jl")
    using Base.Threads
    using JLD2
    include("STMonteCarlo_MultiThread.jl")


#=
    Optimisation of the multi-threaded version would not have been possible without the kind guidence of those on the Julia Forums: https://discourse.julialang.org/t/fast-multi-threaded-array-editing-without-data-races/114863/41
    In particular the assistance of users: mbauman, adienes, Oscar_Smith, Satvik, Salmon, sgaure and foobar_lv2
=#

function SpectraEvaluateMultiThread()


    # ========= Load/Create Files ========== #

        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            SAtotal = f["STotal"];
            TAtotal = f["TTotal"];
            SAtally = f["STally"];
            TAtally = f["TTally"];
            #SMatrix = f["SMatrix"];
            #TMatrix = f["TMatrix"];
            p3Max = f["p3Max"];
            t3MinMax = f["t3MinMax"];
            close(f)
        else
            SAtotal = zeros(Float32,(nump3+1),numt3,nump1,numt1,nump2,numt2); 
            TAtotal = zeros(Float32,nump1,numt1,nump2,numt2);
            SAtally = zeros(UInt32,numt3,nump1,numt1,nump2,numt2);
            TAtally = zeros(UInt32,nump1,numt1,nump2,numt2);
            p3Max = zeros(Float32,numt3,nump1,numt1,nump2,numt2);
            t3MinMax = zeros(Float32,2,(nump3+1),nump1,numt1,nump2,numt2);
        end

    # ====================================== #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:nump1]    

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        # Set up workers
        workers = [STMonteCarloAxi_MultiThread!(SAtotal,TAtotal,SAtally,TAtally,ArrayOfLocks,p3Max,t3MinMax) for i in 1:nThreads]
        
        wait.(workers) # Allow all workers to finish
   
    # ===================================== #

    # ===== Calculate S and T Matricies === #

        #= these operations are run in serial =#

        # preallocate
        SMatrix = zeros(Float32,(nump3+1),numt3,nump1,numt1,nump2,numt2);
        TMatrix = zeros(Float32,nump1,numt1,nump2,numt2);

        # divide element wise by tallys
        for i in axes(SMatrix,1)
            @. @view(SMatrix[i,:,:,:,:,:]) = @view(SAtotal[i,:,:,:,:,:]) / SAtally
        end
        replace!(SMatrix,NaN=>0f0); # remove NaN caused by /0f0
        TMatrix = TAtotal ./ TAtally;
        replace!(TMatrix,NaN=>0f0);

        #testtotal = (dropdims(sum(SAtotal[:,:,60,:,60,:],dims=(3,4)),dims=(3,4)))
        #testtally = Int.(dropdims(sum(SAtally[:,10,:,10,:],dims=(2,3)),dims=(2,3)))
        #testtotal = (dropdims(sum(SMatrix[:,:,40,:,40,:],dims=(3,4)),dims=(3,4)))

        # Angle / Momentum Ranges
        t3val = trange(t3l,t3u,numt3)
        t1val = trange(t1l,t1u,numt1)
        t2val = trange(t2l,t2u,numt2)
        p3val = prange(p3l,p3u,nump3)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)


        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val)    #applies phase space factors for symmetries
        STSymmetry!(SMatrix,TMatrix)                                        #initial states are symmetric -> apply symmetry of interaction to improve MC values
        #PhaseSpaceFactors2!(SMatrix,TMatrix,p3val,t3val,p1val,t1val)    #corrects phase space factors for application in kinetic models
        PhaseSpaceFactors2_3D!(SMatrix,TMatrix,p3val,t3val,p1val,t1val)
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 


    # ===================================== # 

    # ========== Save Arrays ============== #
        
        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal",SAtotal)
        write(f,"TTotal",TAtotal)
        write(f,"STally",SAtally)
        write(f,"TTally",TAtally)
        write(f,"SMatrix",SMatrix)
        write(f,"TMatrix",TMatrix)
        write(f,"p3Max",p3Max)
        write(f,"t3MinMax",t3MinMax)
        write(f,"name1Data",eval(Symbol(name1*"Data")))
        write(f,"name2Data",eval(Symbol(name2*"Data")))
        write(f,"name3Data",eval(Symbol(name3*"Data")))
        write(f,"name4Data",eval(Symbol(name4*"Data")))
        close(f)

    # ===================================== #

    return nothing

end 