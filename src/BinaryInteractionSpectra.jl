module BinaryInteractionSpectra

export SpectraEvaluateSerial, SpectraEvaluateMultiThread, fload, fread, fclose

    using JLD2
    using Base.Threads
    #using AccurateArithmetic
    using KahanSummation

    # include common files
        include("Common/MyPhysicalConstants.jl")
        include("Common/ParticleData.jl")
        include("Common\\Init.jl")
        include("Common\\DifferentialCrossSectionFunctions.jl")
        include("Common\\Momentum3Values.jl")
        include("Common\\RandomPointMomentum.jl")
        include("Common\\RandomPointSphere.jl")
        include("Common/MandelstramChecks.jl")
        include("Common\\STValue.jl")
        include("Common/UsefulGridValueFunctions.jl")
        include("Common/PhaseSpaceFactors.jl")
        include("Common/Location.jl")

    # include serial methods
        include("Serial\\STIntegration_Serial.jl")

    #include parallel methods
        include("MultiThread\\STIntegration_MultiThread.jl")

    function fload()
        
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
            close(f)
        else
            error("no file")
        end

        return (SAtot, TAtot, AStal, ATtal, SMatrix, TMatrix, p3Max, t3MinMax);
        #run (Stot, Ttot, Stal, Ttal, SMatrix, TMatrix, p3Max, t3MinMax) = fload();

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

