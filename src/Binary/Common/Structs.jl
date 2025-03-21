struct BinaryParameters

    name1::String
    name2::String
    name3::String
    name4::String

    mu1::Float64
    mu2::Float64
    mu3::Float64
    mu4::Float64

    p1_low::Float64
    p1_up::Float64
    p1_grid::String
    p1_num::Int64

    u1_grid::String
    u1_num::Int64

    phi1_grid::String
    phi1_num::Int64

    p2_low::Float64
    p2_up::Float64
    p2_grid::String
    p2_num::Int64

    u2_grid::String
    u2_num::Int64

    phi2_grid::String
    phi2_num::Int64

    p3_low::Float64
    p3_up::Float64
    p3_grid::String
    p3_num::Int64

    u3_grid::String
    u3_num::Int64

    phi3_grid::String
    phi3_num::Int64

    p4_low::Float64
    p4_up::Float64
    p4_grid::String
    p4_num::Int64

    u4_grid::String
    u4_num::Int64

    phi4_grid::String
    phi4_num::Int64

end

mutable struct BinaryUserInput <: Function

    Parameters::BinaryParameters
    numTiter::Int64
    numSiter::Int64
    fileLocation::String
    fileName::String
    numThreads::Int64
    MinMax::Bool

    function BinaryUserInput(Parameters::BinaryParameters,fileLoc::String,numT::Int64,numS::Int64;numThreads=1,MinMax=false)
        self = new()

        P = Parameters

        self.Parameters = P
        self.numTiter = numT
        self.numSiter = numS
        self.fileLocation = fileLoc
        self.numThreads = numThreads
        self.MinMax = MinMax

        self.fileName = P.name1*P.name2*P.name3*P.name4
        self.fileName *= "#"*string(P.p1_low)*"-"*string(P.p1_up)*P.p1_grid*string(P.p1_num)
        self.fileName *= "#"*P.u1_grid*string(P.u1_num)
        self.fileName *= "#"*P.phi1_grid*string(P.phi1_num)

        self.fileName *= "#"*string(P.p2_low)*"-"*string(P.p2_up)*P.p2_grid*string(P.p2_num)
        self.fileName *= "#"*P.u2_grid*string(P.u2_num)
        self.fileName *= "#"*P.phi2_grid*string(P.phi2_num)

        self.fileName *= "#"*string(P.p3_low)*"-"*string(P.p3_up)*P.p3_grid*string(P.p3_num)
        self.fileName *= "#"*P.u3_grid*string(P.u3_num)
        self.fileName *= "#"*P.phi3_grid*string(P.phi3_num)

        self.fileName *= "#"*string(P.p4_low)*"-"*string(P.p4_up)*P.p4_grid*string(P.p4_num)
        self.fileName *= "#"*P.u4_grid*string(P.u4_num)
        self.fileName *= "#"*P.phi4_grid*string(P.phi4_num)
        
        self.fileName *= ".jld2";

        return self

    end

end

mutable struct ScatteringArrays <: Function

    SAtotal3::Array{Float64,9}
    SAtotal4::Array{Float64,9}

    SAtally3::Array{Float64,8}
    SAtally4::Array{Float64,8}

    SMatrix3::Array{Float64,9}
    SMatrix4::Array{Float64,9}

    TAtotal::Array{Float64,6}
    TAtally::Array{Float64,6}
    TMatrix1::Array{Float64,6}
    TMatrix2::Array{Float64,6}

    p3Max::Array{Float64,8}
    u3MinMax::Array{Float64,8}
    p4Max::Array{Float64,8}
    u4MinMax::Array{Float64,8}

    function ScatteringArrays(UserInput::BinaryUserInput)
        self = new()

        filePath = UserInput.fileLocation*"\\"*UserInput.fileName
        fileExist = isfile(filePath)
        MinMax = UserInput.MinMax

        if fileExist
            f = jldopen(filePath,"r+");
            self.SAtotal3 = f["STotal3"];
            self.SAtally3 = f["STally3"];
            self.SMatrix3 = f["SMatrix3"];
            self.SAtotal4 = f["STotal4"];
            self.SAtally4 = f["STally4"];
            self.SMatrix4 = f["SMatrix4"];
            self.TAtotal = f["TTotal"];
            self.TAtally = f["TTally"];
            self.TMatrix1 = f["TMatrix1"];
            self.TMatrix2 = f["TMatrix2"];
            if MinMax
                self.p3Max = f["p3Max"];
                self.u3MinMax = f["u3MinMax"];
                self.p4Max = f["p4Max"];
                self.u4MinMax = f["u4MinMax"];
            end
            close(f)
        else
            P = UserInput.Parameters

            p1_num = P.p1_num
            u1_num = P.u1_num
            phi1_num = P.phi1_num

            p2_num = P.p2_num
            u2_num = P.u2_num
            phi2_num = P.phi2_num

            p3_num = P.p3_num
            u3_num = P.u3_num
            phi3_num = P.phi3_num

            p4_num = P.p4_num
            u4_num = P.u4_num
            phi4_num = P.phi4_num

            self.SAtotal3 = zeros(Float64,(p3_num+1),u3_num,phi3_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num); 
            self.SAtotal4 = zeros(Float64,(p4_num+1),u4_num,phi4_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num); 
            self.TAtotal = zeros(Float64,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            self.SAtally3 = zeros(UInt32,u3_num,phi3_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            self.SAtally4 = zeros(UInt32,u4_num,phi4_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            self.TAtally = zeros(UInt32,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            self.SMatrix3 = zeros(Float64,(p3_num+1),u3_num,phi3_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            self.TMatrix1 = zeros(Float64,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            self.TMatrix2 = zeros(Float64,p2_num,u2_num,phi2_num,p1_num,u1_num,phi1_num);
            self.SMatrix4 = zeros(Float64,(p4_num+1),u4_num,phi4_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
            if MinMax
                self.p3Max = zeros(Float64,u3_num,phi3_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
                self.u3MinMax = zeros(Float64,2,(p3_num+1),phi3_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
                self.p4Max = zeros(Float64,u4_num,phi4_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
                self.u4MinMax = zeros(Float64,2,(p4_num+1),phi4_num,p1_num,u1_num,phi1_num,p2_num,u2_num,phi2_num);
                fill!(@view(self.u3MinMax[1,:,:,:,:,:,:,:,:]),1e0);
                fill!(@view(self.u3MinMax[2,:,:,:,:,:,:,:,:]),-1e0);
                fill!(@view(self.u4MinMax[1,:,:,:,:,:,:,:,:]),1e0);
                fill!(@view(self.u4MinMax[2,:,:,:,:,:,:,:,:]),-1e0);
            end
        end

        return self

    end

end