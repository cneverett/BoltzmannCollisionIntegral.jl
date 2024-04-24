module BinaryInteractionSpectra

export SpectraEvaluateSerial, SpectraEvaluateMultiThread
include("Serial\\STIntegration_Serial.jl")
include("MultiThread\\STIntegration_MultiThread.jl")

end
