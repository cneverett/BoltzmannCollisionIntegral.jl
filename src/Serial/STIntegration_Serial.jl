#= Script for running the ST integration and returning data arrays =#

    using JLD2
    include("STMonteCarlo_Serial.jl")
    #include("..\\Common\\UsefulGridValueFunctions.jl")
    #include("..\\Common\\PhaseSpaceFactors.jl")
    
function SpectraEvaluateSerial()

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

    # ========= Pre-Allocate Arrays ======== #

        # pre-allocate arrays for momentum
        p3v = zeros(Float32,3);
        p3pv = zeros(Float32,3); 
        p1v = zeros(Float32,3);
        p2v = zeros(Float32,3);

        # pre-allocate arrays for ST values
        #ST = zeros(Float32,3); # [S,Sp,T]

    # ===================================== #

    # ===== Run MonteCarlo Integration ==== #

        STMonteCarloAxi_Serial!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax)

    # ===================================== #

    # ===== Calculate S and T Matricies === #

        # preallocate
        SMatrixOldSum = dropdims(sum(SMatrix,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix,0f0);
        TMatrixOldSum = dropdims(sum(TMatrix,dims=(3,4)),dims=(3,4));
        fill!(TMatrix,0f0);


        # divide element wise by tallys
        @inbounds for i in axes(SMatrix,1)
            @. @view(SMatrix[i,:,:,:,:,:]) = @view(SAtotal[i,:,:,:,:,:]) / SAtally
        end
        replace!(SMatrix,NaN=>0f0) # remove NaN caused by /0f0
        TMatrix = TAtotal ./ TAtally
        replace!(TMatrix,NaN=>0f0)

        # Angle / Momentum Ranges
        t3val = trange(t3l,t3u,numt3) # bounds of numt3 blocks
        t1val = trange(t1l,t1u,numt1)
        t2val = trange(t2l,t2u,numt2)
        p3val = prange(p3l,p3u,nump3)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val)    #applies phase space factors for symmetries
        STSymmetry!(SMatrix,TMatrix)                                        #initial states are symmetric -> apply symmetry of interaction to improve MC values
        PhaseSpaceFactors2!(SMatrix,TMatrix,p3val,t3val,p1val,t1val)    #corrects phase space factors for application in kinetic models
                                 
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

end #function