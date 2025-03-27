#= Script for running the ST integration and returning data arrays =#

"""
    SpectraEvaluateSerial(userInputSerial)

Function to run the Monte Carlo integration of the S and T arrays in a serial environment. The function will run the Monte Carlo integration in serial and then calculate the S and T matrices and save the results to a file.
"""
function SpectraEvaluateSerial(userInputSerial::BinaryUserInput)

    # ========= Load user Parameters ======= #

        Parameters = userInputSerial.Parameters
        numTiter = userInputSerial.numTiter
        numSiter = userInputSerial.numSiter
        fileLocation = userInputSerial.fileLocation
        fileName = userInputSerial.fileName
        MinMax = userInputSerial.MinMax

        name1 = Parameters.name1
        name2 = Parameters.name2
        name3 = Parameters.name3
        name4 = Parameters.name4
        mu1 = Parameters.mu1
        mu2 = Parameters.mu2
        mu3 = Parameters.mu3
        mu4 = Parameters.mu4

        p1_low = Parameters.p1_low
        p1_up = Parameters.p1_up
        p1_grid = Parameters.p1_grid
        p1_num = Parameters.p1_num

        u1_grid = Parameters.u1_grid
        u1_num = Parameters.u1_num

        h1_grid = Parameters.h1_grid
        h1_num = Parameters.h1_num

        p2_low = Parameters.p2_low
        p2_up = Parameters.p2_up
        p2_grid = Parameters.p2_grid
        p2_num = Parameters.p2_num

        u2_grid = Parameters.u2_grid
        u2_num = Parameters.u2_num

        h2_grid = Parameters.h2_grid
        h2_num = Parameters.h2_num

        p3_low = Parameters.p3_low
        p3_up = Parameters.p3_up
        p3_grid = Parameters.p3_grid
        p3_num = Parameters.p3_num

        u3_grid = Parameters.u3_grid
        u3_num = Parameters.u3_num

        h3_grid = Parameters.h3_grid
        h3_num = Parameters.h3_num

        p4_low = Parameters.p4_low
        p4_up = Parameters.p4_up
        p4_grid = Parameters.p4_grid
        p4_num = Parameters.p4_num

        u4_grid = Parameters.u4_grid
        u4_num = Parameters.u4_num

        h4_grid = Parameters.h4_grid
        h4_num = Parameters.h4_num

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

        println("Loading/Creating Files")

        filePath = fileLocation*"\\"*fileName

        Arrays = ScatteringArrays(userInputSerial)

    # ====================================== #

    # ======= Define Cross Section Functions Based on Particle Selections ======== #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma::Function = getfield(BoltzmannCollisionIntegral,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt::Function = getfield(BoltzmannCollisionIntegral,name_dsigmadt)

    # ============================================================================ #

    # ===== Run MonteCarlo Integration ==== #

        println("Running Monte Carlo Integration")

        STMonteCarloAxi_Serial!(Arrays,sigma,dsigmadt,userInputSerial)

    # ===================================== #

    # ===== Calculate S and T Matrices === #

        println("Building S and T Matrices")

        SMatrix3 = Arrays.SMatrix3
        SAtally3 = Arrays.SAtally3
        SAtotal3 = Arrays.SAtotal3

        SMatrix4 = Arrays.SMatrix4
        SAtally4 = Arrays.SAtally4
        SAtotal4 = Arrays.SAtotal4

        TMatrix1 = Arrays.TMatrix1
        TAtally = Arrays.TAtally
        TAtotal = Arrays.TAtotal
        TMatrix2 = Arrays.TMatrix2

        if MinMax
            p3Max = Arrays.p3Max
            u3MinMax = Arrays.u3MinMax
            p4Max = Arrays.p4Max
            u4MinMax = Arrays.u4MinMax
        end

        # preallocate
        SMatrixOldSum3 = dropdims(sum(SMatrix3,dims=(4,5,6,7,8,9)),dims=(4,5,6,7,8,9));
        fill!(SMatrix3,0e0);
        SMatrixOldSum4 = dropdims(sum(SMatrix4,dims=(4,5,6,7,8,9)),dims=(4,5,6,7,8,9));
        fill!(SMatrix4,0e0);
        TMatrixOldSum = dropdims(sum(TMatrix1,dims=(4,5,6)),dims=(4,5,6));
        fill!(TMatrix1,0e0);

        # divide element wise by tallys
        if Indistinguishable_34 == true
            @. SAtally3 = SAtally3 + SAtally4
            @. SAtotal3 = SAtotal3 + SAtotal4
            for i in axes(SMatrix3,1)
                @. @view(SMatrix3[i,:,:,:,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:,:,:,:]) / SAtally3
            end
            replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            if MinMax
                @. p3Max = max(p3Max,p4Max)
                @view(u3MinMax[1,:,:,:,:,:,:,:,:]) .= min.(u3MinMax[1,:,:,:,:,:,:,:,:],u4MinMax[1,:,:,:,:,:,:,:,:])
                @view(u3MinMax[2,:,:,:,:,:,:,:,:]) .= max.(u3MinMax[2,:,:,:,:,:,:,:,:],u4MinMax[2,:,:,:,:,:,:,:,:])
                # reset arrays to avoid overcounting when multiple runs are made
                fill!(SAtally4,UInt32(0))
                fill!(SAtotal4,0e0)
                fill!(p4Max,0e0)
                fill!(@view(u4MinMax[1,:,:,:,:,:,:,:,:]),1e0);
                fill!(@view(u4MinMax[2,:,:,:,:,:,:,:,:]),-1e0);
            end
        elseif mu3 == mu4 # system symmetric in 34 interchange
            @. SAtally3 = SAtally3 + SAtally4
            @. SAtally4 = SAtally3
            @. SAtotal3 = SAtotal3 + SAtotal4
            @. SAtotal4 = SAtotal3
            for i in axes(SMatrix3,1)
                @. @view(SMatrix3[i,:,:,:,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:,:,:,:]) / SAtally3
            end
            replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            @. SMatrix4 = SMatrix3
            if MinMax
                @. p3Max = max(p3Max,p4Max)
                @. p4Max = p3Max
                @view(u3MinMax[1,:,:,:,:,:,:,:,:]) .= min.(u3MinMax[1,:,:,:,:,:,:,:,:],u4MinMax[1,:,:,:,:,:,:,:,:])
                @view(u3MinMax[2,:,:,:,:,:,:,:,:]) .= max.(u3MinMax[2,:,:,:,:,:,:,:,:],u4MinMax[2,:,:,:,:,:,:,:,:])
                @. u4MinMax = u3MinMax
            end
        else
            for i in axes(SMatrix3,1)
                @. @view(SMatrix3[i,:,:,:,:,:,:,:,:]) = @view(SAtotal3[i,:,:,:,:,:,:,:,:]) / SAtally3
            end
            replace!(SMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            for i in axes(SMatrix4,1)
                @. @view(SMatrix4[i,:,:,:,:,:,:,:,:]) = @view(SAtotal4[i,:,:,:,:,:,:,:,:]) / SAtally4
            end
            replace!(SMatrix4,NaN=>0e0); # remove NaN caused by /0e0
        end
        TMatrix1 = TAtotal ./ TAtally;
        replace!(TMatrix1,NaN=>0e0);

        # Angle / Momentum Ranges
        u3val = bounds(u_low,u_up,u3_num,u3_grid)
        u4val = bounds(u_low,u_up,u4_num,u4_grid)
        u1val = bounds(u_low,u_up,u1_num,u1_grid)
        u2val = bounds(u_low,u_up,u2_num,u2_grid)
        p3val = bounds(p3_low,p3_up,p3_num,p3_grid)
        p4val = bounds(p4_low,p4_up,p4_num,p4_grid)
        p1val = bounds(p1_low,p1_up,p1_num,p1_grid)
        p2val = bounds(p2_low,p2_up,p2_num,p2_grid)
        h3val = bounds(h_low,h_up,h3_num,h3_grid) .* pi
        h4val = bounds(h_low,h_up,h4_num,h4_grid) .* pi
        h1val = bounds(h_low,h_up,h1_num,h1_grid) .* pi
        h2val = bounds(h_low,h_up,h2_num,h2_grid) .* pi

        println("Applying Momentum Space Factors")

        # Momentum space volume elements and symmetries
        PhaseSpaceFactors1!(SMatrix3,SMatrix4,TMatrix1,u3val,h3val,u4val,h4val,p1val,u1val,h1val,p2val,u2val,h2val,Indistinguishable_12)      # applies phase space factors for symmetries                  
        STSymmetry!(SMatrix3,SMatrix4,TMatrix1,mu1,mu2)   # initial states are symmetric -> apply symmetry of interaction to improve MC values
        if Indistinguishable_12 == false
            perm = [3,4,1,2]
            TMatrix2 = permutedims(TMatrix1,perm)
        end
        PhaseSpaceFactors2!(SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3val,u3val,h3val,p4val,u4val,h4val,p1val,u1val,h1val,p2val,u2val,h2val)            # corrects phase space factors for application in kinetic models
                                            
        # correction to better conserve particle number and account for statistical noise of MC method
        #SCorrection2!(SMatrix,TMatrix) 

        # output a measure of convergence, i.e. new-old/old
        SMatrixSum3 = dropdims(sum(SMatrix3,dims=(4,5,6,7,8,9)),dims=(4,5,6,7,8,9));
        SConverge3 = (SMatrixSum3 .- SMatrixOldSum3)./SMatrixOldSum3
        SMatrixSum4 = dropdims(sum(SMatrix4,dims=(4,5,6,7,8,9)),dims=(4,5,6,7,8,9));
        SConverge4 = (SMatrixSum4 .- SMatrixOldSum4)./SMatrixOldSum4
        TMatrixSum = dropdims(sum(TMatrix1,dims=(4,5,6)),dims=(4,5,6));
        TConverge = (TMatrixSum .- TMatrixOldSum)./TMatrixOldSum

    # ===================================== # 

    # ========== Save Arrays ============== #

        #OutputParameters = Parameters

        println("Saving Arrays")

        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"STotal3",SAtotal3)
        write(f,"STally3",SAtally3)
        write(f,"SMatrix3",SMatrix3)
        write(f,"SConverge3",SConverge3)

        write(f,"STotal4",SAtotal4)
        write(f,"STally4",SAtally4)
        write(f,"SMatrix4",SMatrix4)
        write(f,"SConverge4",SConverge4)

        write(f,"TTotal",TAtotal)
        write(f,"TTally",TAtally)
        write(f,"TMatrix1",TMatrix1)
        write(f,"TMatrix2",TMatrix2)

        write(f,"TConverge",TConverge)

        if MinMax
            write(f,"p3Max",p3Max)
            write(f,"u3MinMax",u3MinMax)
            write(f,"p4Max",p4Max)
            write(f,"u4MinMax",u4MinMax)
        end

        write(f,"Parameters",Parameters)
        close(f)

    # ===================================== #

        return nothing

end #function