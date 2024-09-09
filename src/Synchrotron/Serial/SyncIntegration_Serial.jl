"""
    SyncEvaluateSerial!(userSyncInputSerial)

Function to run the Monte Carlo integration of the S array in a serial enviroment. 
"""
function SynEvaluateSerial!(userSyncInputSerial)

    # ======= Load User Parameters ======= #

        (name1,name2,p1l,p1u,nump1,numt1,p2l,p2u,nump2,numt2,numSiter,fileLocation,fileName,BMag) = userSyncInputSerial;

    # ==================================== #

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
            SAtotal = zeros(Float64,nump1,numt1,nump2,numt2);  
            SAtally = zeros(UInt32,nump1,numt1,nump2,numt2);
            SMatrix = zeros(Float64,nump1,numt1,nump2,numt2);
            #pMax = zeros(Float64,nump1,numt1,nump2,numt2);
            #tMinMax = zeros(Float64,2,nump1,numt1,nump2,numt2);
            #fill!(@view(tMinMax[1,:,:,:,:]),1e0);
        end

    # ================================= #

    # ===== Set Particle (normalised) Masses) and Parameters ====== #

        mu2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
        z2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("z"*name2))

        Parameters = (mu2,z2,BMag,p1l,p1u,nump1,p2l,p2u,nump2,numt1,numt2)

    # ============================================================= #


    # ===== Run MonteCarlo Integration ==== #

        STMonteCarloAxi_Serial!(SAtotal,SAtally,#=pMax,tMinMax,=#Parameters,numSiter)

    # ===================================== #


    # ===== Calculate S and T Matricies === #

        # preallocate
        SMatrixOldSum = dropdims(sum(SMatrix,dims=(3,4,5,6)),dims=(3,4,5,6));
        fill!(SMatrix,0e0);
        
        # divide element wise by tallys
        for i in axes(SMatrix,1)
            @. @view(SMatrix[i,:,:,:,:,:]) = @view(SAtotal[i,:,:,:,:,:]) / SAtally
        end
        replace!(SMatrix,NaN=>0e0); # remove NaN caused by /0e0

        # Angle / Momentum Ranges
        t1val = trange(numt1)
        t2val = trange(numt2)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactorsSync1!(SMatrix,p1val,t1val,p2val,t2val)      # applies phase space factors for symmetries                  
        SyncSymmetry!(SMatrix)   # initial states are symmetric -> apply symmetry of interaction to improve MC values
        PhaseSpaceFactorsSync2!(SMatrix,p1val,t1val,p2val,t2val)            # corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

        # output a measure of convergence, i.e. new-old/old
        SMatrixSum = dropdims(sum(SMatrix,dims=(3,4,5,6)),dims=(3,4,5,6));
        SConverge = (SMatrixSum .- SMatrixOldSum)./SMatrixOldSum

    # ===================================== # 

    # ========== Save Arrays ============== #

        OutputParameters = (name1,name2,p1l,p1u,nump1,p2l,p2u,nump2,numt1,numt2)
            
        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal",SAtotal)
        write(f,"STally",SAtally)
        write(f,"SMatrix",SMatrix)
        #write(f,"pMax",pMax)
        #write(f,"tMinMax",tMinMax)
        write(f,"SConverge",SConverge)
        write(f,"Parameters",OutputParameters)
        close(f)

    # ===================================== #

end # function 