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
    #   include("STMonteCarlo_MultiThread.jl")

#=
    Optimisation of the multi-threaded version would not have been possible without the kind guidence of those on the Julia Forums: https://discourse.julialang.org/t/fast-multi-threaded-array-editing-without-data-races/114863/41
    In particular the assistance of users: mbauman, adienes, Oscar_Smith, Satvik, Salmon, sgaure and foobar_lv2
=#

"""
    STMonteCarloAxi_MultiThread!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax})

Intput:

    - Domain Boundaries (defined as CONST in Init.jl)
        - p bounds and divisions for species 1,3,4
        - theta divisions for species 1,3,4 ( bounds not nessessarlity needed as assumed [0,1] )
        - phi divisions for species 1,3,4 ( bounds not nessessarlity needed as assumed [0,2] )
    - Particle Masses (defined as CONST in Init.jl)
        - for species 1,2,3,4
    - Array of stored integration totals and tallys 
        - total is cumulative sum of reaction rate in that domain
        - tally is cumalitive total of points that have been sampled in that doimain
        - S Array will have dimensions ((nump3+1) x numt3 x nump1 x numt1 x nump2 x numt2) for axisymmetric
            - extra entry for p3 is for overflow momenta i.e. array acts like [p3 i, p3 i+1, p3 i+2 .... p3 nump3, overflow]
        - T Array will have dimensions (nump1 x numt1 x nump2 x numt2) for axisymmetric
    - numTiter and numSiter (defined in Init.jl) as the number of T and S integrations to perform.
    - nThreads as the number of threads to run the integration on

Calculation:

    - Set up workers to perform the integration on multiple threads
    - Random Sample points in each of these domains
        - RandomPointSphere for theta and phi (for species 1,2,3)
        - RandomPointMomentum for p ( species 1,2 only )
    - Take random points (t3,h1,p1,p2,t1,t2,h3,h4) and calculate valid p3 point/points 
    - Find position in S and T arrays and allocated tallies and totals accordingly (using locks to ensure single thread access to arrays). 

Output:

    - Edited arrays of stored integration totals and tallys
    - One array (S array) gives rate of reaction to particular state/particles 1(2) from state 34 i.e. rate of emission of 1 from reaction 34->1(2)
    - One array (T array) gives rate of reaction from state/particles 34 to any state 12 i.e. rate of absorption of 34 in reaction 34->12

"""
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
            SMatrix = f["SMatrix"];
            TMatrix = f["TMatrix"];
            p3Max = f["p3Max"];
            t3MinMax = f["t3MinMax"];
            close(f)
        else
            SAtotal = zeros(Float32,(nump3+1),numt3,nump1,numt1,nump2,numt2); 
            TAtotal = zeros(Float32,nump1,numt1,nump2,numt2);
            SAtally = zeros(UInt32,numt3,nump1,numt1,nump2,numt2);
            TAtally = zeros(UInt32,nump1,numt1,nump2,numt2);
            SMatrix = zeros(Float32,(nump3+1),numt3,nump1,numt1,nump2,numt2);
            TMatrix = zeros(Float32,nump1,numt1,nump2,numt2);
            p3Max = zeros(Float32,numt3,nump1,numt1,nump2,numt2);
            t3MinMax = zeros(Float32,2,(nump3+1),nump1,numt1,nump2,numt2);
        end

    # ====================================== #

    # ======== Set up Array of Locks ====== #

        ArrayOfLocks = [Threads.SpinLock() for _ in 1:nump1]    

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        # Set up workers
        workers = [STMonteCarloAxi_MultiThread!(SAtotal,TAtotal,SAtally,TAtally,ArrayOfLocks,p3Max,t3MinMax) for _ in 1:nThreads]
        
        wait.(workers) # Allow all workers to finish
   
    # ===================================== #

    # ===== Calculate S and T Matricies === #

        #= these operations are run in serial =#

        # preallocate
        SMatrixOldSum = dropdims(sum(SMatrix,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix,0f0);
        TMatrixOldSum = dropdims(sum(TMatrix,dims=(3,4)),dims=(3,4));
        fill!(TMatrix,0f0);

        # divide element wise by tallys
        for i in axes(SMatrix,1)
            @. @view(SMatrix[i,:,:,:,:,:]) = @view(SAtotal[i,:,:,:,:,:]) / SAtally
        end
        replace!(SMatrix,NaN=>0f0); # remove NaN caused by /0f0
        TMatrix = TAtotal ./ TAtally;
        replace!(TMatrix,NaN=>0f0);

        # Angle / Momentum Ranges
        t3val = trange(t3l,t3u,numt3)
        t1val = trange(t1l,t1u,numt1)
        t2val = trange(t2l,t2u,numt2)
        p3val = prange(p3l,p3u,nump3)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val)      # applies phase space factors for symmetries                  
        STSymmetry!(SMatrix,TMatrix)                                            # initial states are symmetric -> apply symmetry of interaction to improve MC values
        PhaseSpaceFactors2!(SMatrix,TMatrix,p3val,t3val,p1val,t1val)            # corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

        # output a measure of convergence, i.e. new-old/old
        SMatrixSum = dropdims(sum(SMatrix,dims=(3,4,5,6)),dims=(3,4,5,6));
        SConverge = (SMatrixSum .- SMatrixOldSum)./SMatrixOldSum
        TMatrixSum = dropdims(sum(TMatrix,dims=(3,4)),dims=(3,4));
        TConverge = (TMatrixSum .- TMatrixOldSum)./TMatrixOldSum


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
        write(f,"SConverge",SConverge)
        write(f,"TConverge",TConverge)
        write(f,"name1Data",eval(Symbol(name1*"Data")))
        write(f,"name2Data",eval(Symbol(name2*"Data")))
        write(f,"name3Data",eval(Symbol(name3*"Data")))
        write(f,"name4Data",eval(Symbol(name4*"Data")))
        close(f)

    # ===================================== #

    return nothing

end # function