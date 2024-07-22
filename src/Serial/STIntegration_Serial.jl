#= Script for running the ST integration and returning data arrays =#

"""
    SpectraEvaluateSerial(userInputSerial)

Function to run the Monte Carlo integration of the S and T arrays in a serial environment. The function will run the Monte Carlo integration in serial and then calculate the S and T matricies and save the results to a file.
"""
function SpectraEvaluateSerial(userInputSerial::Tuple{String,String,String,String,Float32,Float32,Int64,Float32,Float32,Int64,Float32,Float32,Int64,Int64,Int64,Int64,Int64,Int64,String,String})

    # ========= Load user Parameters ======= #

        (name1,name2,name3,name4,p3l,p3u,nump3,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt1,numt2,numTiter,numSiter,fileLocation,fileName) = userInputSerial

    # ====================================== #

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

    # ======= Define Cross Section Functions Based on Particle Selections ======== #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma::Function = getfield(BoltzmannCollisionIntegral,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt::Function = getfield(BoltzmannCollisionIntegral,name_dsigmadt)

    # ============================================================================ #

    # ===== Set Particle (normalised) Masses) and Parameters ====== #

        mu1::Float32 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
        mu2::Float32 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
        mu3::Float32 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
        mu4::Float32 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

        Parameters = (mu1,mu2,mu3,mu4,p3l,p3u,nump3,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt1,numt2)

    # ============================================================= #

    # ===== Run MonteCarlo Integration ==== #

        STMonteCarloAxi_Serial!(SAtotal,TAtotal,SAtally,TAtally,p3v,p3pv,p1v,p2v,p3Max,t3MinMax,sigma,dsigmadt,Parameters,numTiter,numSiter)

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
        t3val = trange(numt3) # bounds of numt3 blocks
        t1val = trange(numt1)
        t2val = trange(numt2)
        p3val = prange(p3l,p3u,nump3)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix,TMatrix,t3val,p1val,t1val,p2val,t2val,name1,name2)    #applies phase space factors for symmetries
        STSymmetry!(SMatrix,TMatrix,mu1,mu2)                                        #initial states are symmetric -> apply symmetry of interaction to improve MC values
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

        OutputParameters = (name1,name2,name3,name4,p3l,p3u,nump3,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt1,numt2)
        
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
        write(f,"Parameters",OutputParameters)
        close(f)

    # ===================================== #

        return nothing

end #function