module BoltzmannCollisionIntegral

export SpectraEvaluateSerial, SpectraEvaluateMultiThread, fload_All, DoesConserve, fload_Matrix, fload_Matrix_ISO
export SyncEvaluateSerial, SyncEvaluateMultiThread, fload_All_Sync, fload_Matrix_Sync, fload_Matrix_SyncISO

    using JLD2
    using Base.Threads
    #using StaticArrays
    using BenchmarkTools
    #using Documenter
    using Bessels
    using ProgressMeter

    # include Common files
        include("Common/Constants.jl")
        include("Common/RandomPoints.jl")
        include("Common/Location.jl")

    # include Binary files
        include("Binary/Common/Structs.jl")
        include("Binary/Common/DifferentialCrossSectionFunctions.jl")
        include("Binary/Common/Momentum3Values.jl")
        include("Binary/Common/MandelstramChecks.jl")
        include("Binary/Common/STValue.jl")
        include("Binary/Common/UsefulGridValueFunctions.jl")
        include("Binary/Common/PhaseSpaceFactors.jl")
        include("Binary/Common/Sampling.jl")
        # include serial methods
        include("Binary/Serial/STIntegration_Serial.jl")
        include("Binary/Serial/STMonteCarlo_Serial.jl")
        #include parallel methods
        include("Binary/MultiThread/STIntegration_MultiThread.jl")
        include("Binary/MultiThread/STMonteCarlo_MultiThread.jl")       
        # include data reading functions for export
        include("Binary/Common/DataReading.jl")
        
    # include Synchrotron functions
        include("Emission/Common/SynchrotronKernel.jl")
        include("Emission/Common/SyncPhaseSpaceFactors.jl")
        # include serial methods
        include("Emission/Serial/SyncMonteCarlo_Serial.jl")
        include("Emission/Serial/SyncIntegration_Serial.jl")
        # include parallel methods
        include("Emission/MultiThread/SyncMonteCarlo_MultiThread.jl")
        include("Emission/MultiThread/SyncIntegration_MultiThread.jl")
        # include data reading functions for export
        include("Emission/Common/SyncDataReading.jl")

end

