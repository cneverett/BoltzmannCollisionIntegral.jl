module BinaryInteractionSpectra

export SpectraEvaluateSerial, SpectraEvaluateMultiThread

    using JLD2

    # include common files
        include("Common/MyPhysicalConstants.jl")
        include("Common/ParticleData.jl")
        include("Common\\Init.jl")
        include("Common\\DifferentialCrossSectionFunctions.jl")
        include("Common\\Momentum3Values.jl")
        include("Common\\RandomPointMomentum.jl")
        include("Common\\RandomPointSphere.jl")
        include("Common\\STValue.jl")
        include("Common/UsefulGridValueFunctions.jl")
        include("Common/PhaseSpaceFactors.jl")
        include("Common/location.jl")

    # include serial methods
        include("Serial\\STIntegration_Serial.jl")

    #include parallel methods
        include("MultiThread\\STIntegration_MultiThread.jl")

end

