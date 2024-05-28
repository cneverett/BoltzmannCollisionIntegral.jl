#= Script for running the ST integration and returning data arrays =#

    include("STMonteCarlo_MultiThread.jl")
    #include("..\\Common\\UsefulGridValueFunctions.jl")
    #include("..\\Common\\PhaseSpaceFactors.jl")
    using JLD2

function SpectraEvaluateMultiThread()

    # ========= Load/Create Files ========== #

        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            SAtot = f["STotal"];
            TAtot = f["TTotal"];
            AStal = f["STally"];
            ATtal = f["TTally"];
            #SMatrix = f["SMatrix"];
            #TMatrix = f["TMatrix"];
            close(f)
        else
            SAtot = zeros(Float32,(nump3+2),numt3,nump1,numt1,nump2,numt2); 
            TAtot = zeros(Float32,nump1,numt1,nump2,numt2);
            AStal = zeros(UInt32,(nump3+2),numt3,nump1,numt1,nump2,numt2);
            ATtal = zeros(UInt32,nump1,numt1,nump2,numt2);
        end

    # ====================================== #

    # ========= Pre-Allocate Arrays ======== #

        # pre-allocate arrays for momentum
        p3v = zeros(Float32,3,2,nThreads); # two three array vector ((p3,t3,h1),(p3',t3',h1')) second corresponds to mirrored point in angle space
        p1v = zeros(Float32,3,nThreads);
        p2v = zeros(Float32,3,nThreads);

        # pre-allocate arrays for ST values
        ST = zeros(Float32,3,nThreads); # [S,Sp,T]

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        STMonteCarloAxi_MultiThread!(SAtot,TAtot,AStal,ATtal,p3v,p1v,p2v,ST)

    # ===================================== #

    # ===== Calculate S and T Matricies === #

        #= these operations are run in serial =#

        # preallocate
        SMatrix = zeros(Float32,(nump3+2),numt3,nump1,numt1,nump2,numt2);
        TMatrix = zeros(Float32,nump1,numt1,nump2,numt2);

        # divide element wise by tallys
        SMatrix = SAtot ./ AStal;
        replace!(SMatrix,NaN=>0f0); # remove NaN caused by /0f0
        TMatrix = TAtot ./ ATtal;
        replace!(TMatrix,NaN=>0f0);

        # Angle / Momentum Ranges
        t3val = trange(t3l,t3u,numt3); # bounds of numt3 blocks
        t1val = trange(t1l,t1u,numt1);
        t2val = trange(t2l,t2u,numt2);
        #if (log10pspace == true)
            p3val = prange(p3l,p3u,nump3);
            p1val = prange(p1l,p1u,nump1);
            p2val = prange(p2l,p2u,nump2);
        #= elseif (log10pspace == false)
            p3val = range(p3l,p3u,nump3);
            p1val = range(p1l,p1u,nump1);
            p2val = range(p2l,p2u,nump2);
        else
            error("Log10pspace not defined")
        end;
        =#

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix,TMatrix,p3val,t3val,p1val,t1val,p2val,t2val)    #applies phase space factors for symmetries
        STSymmetry!(SMatrix,TMatrix,mu1,mu2)                                        #initial states are symmetric -> apply symmetry of interaction to improve MC values
        PhaseSpaceFactors2!(SMatrix,TMatrix,p3val,t3val,p1val,t1val,p2val,t2val)    #corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 


    # ===================================== # 

    # ========== Save Arrays ============== #
        
        # have to delete data field and recreate cannot just update
        if fileExist    # i.e. not first time
            f = jldopen(filePath,"r+")
            Base.delete!(f,"STotal")
            Base.delete!(f,"TTotal") 
            Base.delete!(f,"STally")
            Base.delete!(f,"TTally")
            Base.delete!(f,"SMatrix")
            Base.delete!(f,"TMatrix")
            write(f,"STotal",SAtot)
            write(f,"TTotal",TAtot)
            write(f,"STally",AStal)
            write(f,"TTally",ATtal)
            write(f,"SMatrix",SMatrix)
            write(f,"TMatrix",TMatrix)
        else    # create file
            f = jldopen(filePath,"w") # creates file
            write(f,"STotal",SAtot)
            write(f,"TTotal",TAtot)
            write(f,"STally",AStal)
            write(f,"TTally",ATtal)
            write(f,"SMatrix",SMatrix)
            write(f,"TMatrix",TMatrix)
        end
        close(f)

        # --------- Saving Integration Parameters ------ #

        if fileExist==false # only on first time
            f = jldopen(filePath,"r+");
            write(f,"name1Data",eval(Symbol(name1*"Data")))
            write(f,"name2Data",eval(Symbol(name2*"Data")))
            write(f,"name3Data",eval(Symbol(name3*"Data")))
            write(f,"name4Data",eval(Symbol(name4*"Data")))
            close(f)
        end

        # ---------------------------------------------- #

    # ===================================== #


end 