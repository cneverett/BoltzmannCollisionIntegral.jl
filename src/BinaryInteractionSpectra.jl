module BinaryInteractionSpectra

export SpectraEvaluateSerial, SpectraEvaluateMultiThread, fload_All, fread, fclose

    using JLD2
    using Base.Threads
    #using StaticArrays
    using BenchmarkTools
    using Documenter

    # include common files
        include("Common/Constants.jl")
        include("Common/DifferentialCrossSectionFunctions.jl")
        include("Common/Momentum3Values.jl")
        include("Common/RandomPoints.jl")
        include("Common/MandelstramChecks.jl")
        include("Common/STValue.jl")
        include("Common/UsefulGridValueFunctions.jl")
        include("Common/PhaseSpaceFactors.jl")
        include("Common/Location.jl")

    # include serial methods
        include("Serial/STIntegration_Serial.jl")
        include("Serial/STMonteCarlo_Serial.jl")

    #include parallel methods
        include("MultiThread/STIntegration_MultiThread.jl")
        include("MultiThread/STMonteCarlo_MultiThread.jl")

    function fload_All(fileLocation::String,fileName::String)
        
        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r+");
            SAtot = f["STotal"];
            TAtot = f["TTotal"];
            AStal = f["STally"];
            ATtal = f["TTally"];
            SMatrix = f["SMatrix"];
            TMatrix = f["TMatrix"];
            p3Max = f["p3Max"];
            t3MinMax = f["t3MinMax"];
            SConv = f["SConverge"];
            TConv = f["TConverge"]
            close(f)
        else
            error("no file")
        end

        return (SAtot, TAtot, AStal, ATtal, SMatrix, TMatrix, p3Max, t3MinMax,SConv,TConv);
        #run (Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName); in REPL

    end

    function fread()
        filePath = fileLocation*"\\"*fileName
        fileExist = isfile(filePath)

        if fileExist
            f = jldopen(filePath,"r");
            return f
        else
            error("no file")
        end
    end

end

