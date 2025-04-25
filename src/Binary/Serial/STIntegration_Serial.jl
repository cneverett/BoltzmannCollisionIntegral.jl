#= Script for running the ST integration and returning data arrays =#

"""
    SpectraEvaluateSerial(userInputSerial)

Function to run the Monte Carlo integration of the S and T arrays in a serial environment. The function will run the Monte Carlo integration in serial and then calculate the S and T matrices and save the results to a file.
"""
function SpectraEvaluateSerial(userInputSerial::Tuple{Tuple{String,String,String,String,Float64,Float64,Float64,Float64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64, Float64,Float64,String,Int64,String,Int64,String,Int64},Int64,Int64,String,String}#=userInputSerial::BinaryUserInput=#)

    # ========= Load user Parameters ======= #

        (Parameters,numTiter,numSiter,fileLocation,fileName) = userInputSerial
        (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

        #=Parameters = userInputSerial.Parameters
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
        h4_num = Parameters.h4_num=#

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

    # ==== Building Monte Carlo Arrays ===== #

        println("Building Monte Carlo Arrays")

        MCArrays = MonteCarloArrays(Parameters)

        GainTally3 = MCArrays.GainTally3
        GainTotal3 = MCArrays.GainTotal3
        GainTally4 = MCArrays.GainTally4
        GainTotal4 = MCArrays.GainTotal4
        GainMatrix3 = MCArrays.GainMatrix3
        GainMatrix4 = MCArrays.GainMatrix4

        LossTally = MCArrays.LossTally
        LossTotal = MCArrays.LossTotal
        LossMatrix1 = MCArrays.LossMatrix1
        LossMatrix2 = MCArrays.LossMatrix2

    # ====================================== #

    # ======= Define Cross Section Functions Based on Particle Selections ======== #

        name_sigma = Symbol("sigma_"*name1*name2*name3*name4)
        sigma::Function = getfield(BoltzmannCollisionIntegral,name_sigma)
        name_dsigmadt = Symbol("dsigmadt_"*name1*name2*name3*name4)
        dsigmadt::Function = getfield(BoltzmannCollisionIntegral,name_dsigmadt)

    # ============================================================================ #

    # ===== Run MonteCarlo Integration ==== #

        println("Running Monte Carlo Integration")

        STMonteCarlo_Serial!(GainTotal3,GainTotal4,LossTotal,GainTally3,GainTally4,LossTally,sigma,dsigmadt,Parameters,numTiter,numSiter)

    # ===================================== #

    # === Update Gain and Loss Matrices === #

        println("Loading and Updating Gain and Loss Matrices")
        
        filePath = fileLocation*"\\"*fileName

        OldMCArrays = OldMonteCarloArrays(Parameters,filePath)
        OldGainMatrix3 = OldMCArrays.GainMatrix3
        OldGainMatrix4 = OldMCArrays.GainMatrix4
        OldGainTally3 = OldMCArrays.GainTally3
        OldGainTally4 = OldMCArrays.GainTally4
        OldLossMatrix1 = OldMCArrays.LossMatrix1
        OldLossMatrix2 = OldMCArrays.LossMatrix2
        OldLossTally = OldMCArrays.LossTally

        # N values are last element of the tally array
        GainTally3_N = @view(GainTally3[end,:,:,:,:,:,:,:,:])
        GainTally4_N = @view(GainTally4[end,:,:,:,:,:,:,:,:])
        # K value are all but last element of the tally array
        GainTally3_K = @view(GainTally3[1:end-1,:,:,:,:,:,:,:,:])
        GainTally4_K = @view(GainTally4[1:end-1,:,:,:,:,:,:,:,:])

        println("Applying Symmetries")

        # Apply Symmetries to the Gain and Loss Totals and Tallies
        GainLossSymmetry!(GainTotal3,GainTotal4,GainTally3,GainTally4,LossTotal,LossTally,mu1,mu2,mu3,mu4)

        println("Calculating New Gain and Loss Matrices")

        # calculate the gain and loss matrices
        if Indistinguishable_34 == true
            # Only need to calculate GainMatrix3
            for i in axes(GainTotal3,1)
                @. @view(GainMatrix3[i,:,:,:,:,:,:,:,:]) = @view(GainTotal3[i,:,:,:,:,:,:,:,:]) / GainTally3_N
            end
            replace!(GainMatrix3,NaN=>0e0); # remove NaN caused by / 0
            if mu3 == mu4 
                @. GainMatrix4 = GainMatrix3
            else
                GainMatrix4 = zeros(similar(GainMatrix3))
            end
        else
            for i in axes(GainTotal3,1)
                @. @view(GainMatrix3[i,:,:,:,:,:,:,:,:]) = @view(GainTotal3[i,:,:,:,:,:,:,:,:]) / GainTally3_N
            end
            replace!(GainMatrix3,NaN=>0e0); # remove NaN caused by /0e0
            for i in axes(GainTotal4,1)
                @. @view(GainMatrix4[i,:,:,:,:,:,:,:,:]) = @view(GainTotal4[i,:,:,:,:,:,:,:,:]) / GainTally4_N
            end
            replace!(GainMatrix4,NaN=>0e0); # remove NaN caused by /0e0
        end
        LossMatrix1 = LossTotal ./ LossTally;
        replace!(LossMatrix1,NaN=>0e0);
        

        println("Applying Momentum Space Factors")

        # Angle / Momentum Ranges
        u3val = bounds(u_low,u_up,u3_num,u3_grid)
        u4val = bounds(u_low,u_up,u4_num,u4_grid)
        h3val = bounds(h_low,h_up,h3_num,h3_grid).*pi
        h4val = bounds(h_low,h_up,h4_num,h4_grid).*pi

        # Momentum space volume elements and symmetries
        MomentumSpaceFactorsNew!(GainMatrix3,GainMatrix4,u3val,h3val,u4val,h4val,Indistinguishable_12)
                                    
        println("Weighting average of New and Old Gain and Loss Matrices")

        WeightedAverageGain!(GainMatrix3,OldGainMatrix3,GainTally3_K,OldGainTally3,GainMatrix4,OldGainMatrix4,GainTally4_K,OldGainTally4)
        
        WeightedAverageLoss!(LossMatrix1,OldLossMatrix1,LossTally,OldLossTally)

        if Indistinguishable_12 == false
            perm = [4,5,6,1,2,3]
            OldLossMatrix2 = permutedims(OldLossMatrix1,perm)
        else
            OldLossMatrix2 = zeros(similar(OldLossMatrix1))
        end

    # ===================================== #
    
    # ============= Error Estimates ======= # 

        println("Calculating Error Estimates")

        ErrorOutput =  DoesConserve2((Parameters,OldGainMatrix3,OldGainMatrix4,OldLossMatrix1,OldLossMatrix2))

    # ===================================== #

    # ========== Save Arrays ============== #

        println("Saving Arrays")

        f = jldopen(filePath,"w") # creates file and overwrites previous file if one existed
        write(f,"GainTally3",OldGainTally3)
        write(f,"GainMatrix3",OldGainMatrix3)

        write(f,"GainTally4",OldGainTally4)
        write(f,"GainMatrix4",OldGainMatrix4)

        write(f,"LossTally",OldLossTally)
        write(f,"LossMatrix1",OldLossMatrix1)
        write(f,"LossMatrix2",OldLossMatrix2)

        write(f,"Parameters",Parameters)
        write(f,"ErrorEstimates",ErrorOutput)
        close(f)

    # ===================================== #

    return nothing

end #function