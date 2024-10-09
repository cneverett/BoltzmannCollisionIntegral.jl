#= Script for running the ST integration and returning data arrays =#

"""
    SpectraEvaluateSerial(userInputSerial)

Function to run the Monte Carlo integration of the S and T arrays in a serial environment. The function will run the Monte Carlo integration in serial and then calculate the S and T matricies and save the results to a file.
"""
function SpectraEvaluateSerial(userInputSerial::Tuple{String,String,String,String,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Float64,Float64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,String,String})

    # ========= Load user Parameters ======= #

        (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2,numTiter,numSiter,fileLocation,fileName) = userInputSerial

    # ====================================== #

    # ===== Are states Distinguishable ===== #

        Indistinguishable_12 = name1 == name2
        Indistinguishable_34 = name3 == name4

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
            t3MinMax = f["t3MinMax"];
            SAtotal4 = f["STotal4"];
            SAtally4 = f["STally4"];
            SMatrix4 = f["SMatrix4"];
            p4Max = f["p4Max"];
            t4MinMax = f["t4MinMax"];
            TAtotal = f["TTotal"];
            TAtally = f["TTally"];
            TMatrix1 = f["TMatrix1"];
            TMatrix2 = f["TMatrix2"];
            close(f)
        else
            SAtotal3 = zeros(Float64,(nump3+1),numt3,nump1,numt1,nump2,numt2); 
            SAtotal4 = zeros(Float64,(nump4+1),numt4,nump1,numt1,nump2,numt2); 
            TAtotal = zeros(Float64,nump1,numt1,nump2,numt2);
            SAtally3 = zeros(UInt32,numt3,nump1,numt1,nump2,numt2);
            SAtally4 = zeros(UInt32,numt4,nump1,numt1,nump2,numt2)
            TAtally = zeros(UInt32,nump1,numt1,nump2,numt2);
            SMatrix3 = zeros(Float64,(nump3+1),numt3,nump1,numt1,nump2,numt2);
            TMatrix1 = zeros(Float64,nump1,numt1,nump2,numt2);
            TMatrix2 = zeros(Float64,nump2,numt2,nump1,numt1);
            p3Max = zeros(Float64,numt3,nump1,numt1,nump2,numt2);
            t3MinMax = zeros(Float64,2,(nump3+1),nump1,numt1,nump2,numt2);
            SMatrix4 = zeros(Float64,(nump4+1),numt4,nump1,numt1,nump2,numt2);
            p4Max = zeros(Float64,numt4,nump1,numt1,nump2,numt2);
            t4MinMax = zeros(Float64,2,(nump4+1),nump1,numt1,nump2,numt2);
            fill!(@view(t3MinMax[1,:,:,:,:,:]),1e0);
            fill!(@view(t3MinMax[2,:,:,:,:,:]),-1e0);
            fill!(@view(t4MinMax[1,:,:,:,:,:]),1e0);
            fill!(@view(t4MinMax[2,:,:,:,:,:]),-1e0);
        end

    # ====================================== #

    # ======= Define Cross Section Functions Based on Particle Selections ======== #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma::Function = getfield(BoltzmannCollisionIntegral,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt::Function = getfield(BoltzmannCollisionIntegral,name_dsigmadt)

    # ============================================================================ #

    # ===== Set Particle (normalised) Masses) and Parameters ====== #

        mu1::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
        mu2::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
        mu3::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
        mu4::Float64 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

        Parameters = (mu1,mu2,mu3,mu4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2)

    # ============================================================= #

    # ===== Run MonteCarlo Integration ==== #

        STMonteCarloAxi_Serial!(SAtotal3,SAtotal4,TAtotal,SAtally3,SAtally4,TAtally,p3Max,p4Max,t3MinMax,t4MinMax,sigma,dsigmadt,Parameters,numTiter,numSiter)

    # ===================================== #

    # ===== Calculate S and T Matricies === #

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
            @view(t3MinMax[1,:,:,:,:,:]) .= min.(t3MinMax[1,:,:,:,:,:],t4MinMax[1,:,:,:,:,:])
            @view(t3MinMax[2,:,:,:,:,:]) .= max.(t3MinMax[2,:,:,:,:,:],t4MinMax[2,:,:,:,:,:])
            # reset arrays to avoid overcounting when multiple runs are made
                fill!(SAtally4,UInt32(0))
                fill!(SAtotal4,0e0)
                fill!(p4Max,0e0)
                fill!(@view(t4MinMax[1,:,:,:,:,:]),1e0);
                fill!(@view(t4MinMax[2,:,:,:,:,:]),-1e0);
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
            @view(t3MinMax[1,:,:,:,:,:]) .= min.(t3MinMax[1,:,:,:,:,:],t4MinMax[1,:,:,:,:,:])
            @view(t3MinMax[2,:,:,:,:,:]) .= max.(t3MinMax[2,:,:,:,:,:],t4MinMax[2,:,:,:,:,:])
            @. t4MinMax = t3MinMax
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
        t3val = trange(numt3)
        t4val = trange(numt4)
        t1val = trange(numt1)
        t2val = trange(numt2)
        p3val = prange(p3l,p3u,nump3)
        p4val = prange(p4l,p4u,nump4)
        p1val = prange(p1l,p1u,nump1)
        p2val = prange(p2l,p2u,nump2)

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix3,SMatrix4,TMatrix1,t3val,t4val,p1val,t1val,p2val,t2val,Indistinguishable_12)      # applies phase space factors for symmetries                  
        STSymmetry!(SMatrix3,SMatrix4,TMatrix1,mu1,mu2)   # initial states are symmetric -> apply symmetry of interaction to improve MC values
        if Indistinguishable_12 == false
            perm = [3,4,1,2]
            TMatrix2 = permutedims(TMatrix1,perm)
        end
        PhaseSpaceFactors2!(SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3val,t3val,p4val,t4val,p1val,t1val,p2val,t2val)            # corrects phase space factors for application in kinetic models
                                            
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

        OutputParameters = (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2)
            
       f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal3",SAtotal3)
        write(f,"STally3",SAtally3)
        write(f,"SMatrix3",SMatrix3)
        write(f,"p3Max",p3Max)
        write(f,"t3MinMax",t3MinMax)
        write(f,"SConverge3",SConverge3)

        write(f,"STotal4",SAtotal4)
        write(f,"STally4",SAtally4)
        write(f,"SMatrix4",SMatrix4)
        write(f,"p4Max",p4Max)
        write(f,"t4MinMax",t4MinMax)
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

end #function