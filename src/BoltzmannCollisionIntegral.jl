module BoltzmannCollisionIntegral

export SpectraEvaluateSerial, SpectraEvaluateMultiThread, fload_All, DoesConserve, fload_Matrix, fload_Matrix_ISO
export SyncEvaluateSerial, SyncEvaluateMultiThread, fload_All_Sync, fload_Matrix_Sync, fload_Matrix_SyncISO

    using JLD2
    using Base.Threads
    #using StaticArrays
    using BenchmarkTools
    #using Documenter
    using Bessels

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

    # include data reading functions for export
        include("Common/DataReading.jl")
        
    # include Synchrotron functions
        include("Synchrotron/Common/SynchrotronKernel.jl")
        include("Synchrotron/Common/SyncPhaseSpaceFactors.jl")
        include("Synchrotron/Serial/SyncMonteCarlo_Serial.jl")
        include("Synchrotron/Serial/SyncIntegration_Serial.jl")
        include("Synchrotron/MultiThread/SyncMonteCarlo_MultiThread.jl")
        include("Synchrotron/MultiThread/SyncIntegration_MultiThread.jl")
        include("Synchrotron/Common/SyncDataReading.jl")

end

