"""
    fload_All(fileLocation,fileName)

Loads all the data stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    (Parameters,SAtot3,SAtot4,TAtot,SAtal3,SAtal4,TAtal,SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3Max,p4Max,u3MinMax,u4MinMax,SConv3,SConv4,TConv) = fload_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `Stot3` : A 6D matrix totalling all the emission spectrum values sampled for 12->34 interaction.
- `Stot4` : A 6D matrix totalling all the emission spectrum values sampled for 12->43 interaction.
- `Ttot` : A 4D matrix totalling all the absorption spectrum values sampled.
- `Stal3` : A 5D matrix of tallies of the number of emission spectrum values sampled for 12->34 interaction.
- `Stal4` : A 5D matrix of tallies of the number of emission spectrum values sampled for 12->43 interaction.
- `Ttal` : A 4D matrix of tallies of the number of absorption spectrum values sampled.
- `SMatrix3` : A 6D matrix of the emission spectrum for 12->34 interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for 12->43 interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for 12->34 interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for 21->34 interaction i.e. by permutation of TMatrix1 and correct application of phase space factors if species 1 != species 2.
- `p3Max` : The maximum value of the momentum space variable p3 sampled for each bin. (Useful for correcting numerical diffusion)
- `u3MinMax` : The minimum and maximum values of the momentum space variable t3 sampled for each bin. (Useful for correcting numerical diffusion)
- `p4Max` : The maximum value of the momentum space variable p4 sampled for each bin. (Useful for correcting numerical diffusion)
- `u4MinMax` : The minimum and maximum values of the momentum space variable t4 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv3` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Parameters` for 12->34 interaction.
- `SConv4` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Parameters` for 12->43 interaction.
- `TConv` : A 4D matrix of the convergence of the absorption spectrum compared to the previous run with given `Parameters`.

"""
function fload_All(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]

        SAtot3 = f["STotal3"];
        SAtal3 = f["STally3"];
        SMatrix3 = f["SMatrix3"];
        SConv3 = f["SConverge3"];
        p3Max = f["p3Max"];
        u3MinMax = f["u3MinMax"];

        SAtot4 = f["STotal4"];
        SAtal4 = f["STally4"];
        SMatrix4 = f["SMatrix4"];
        SConv4 = f["SConverge4"];
        p4Max = f["p4Max"];
        u4MinMax = f["u4MinMax"];

        TAtot = f["TTotal"];
        TAtal = f["TTally"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];
        TConv = f["TConverge"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,SAtot3,SAtot4,TAtot,SAtal3,SAtal4,TAtal,SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3Max,p4Max,u3MinMax,u4MinMax,SConv3,SConv4,TConv);

end


"""
    fload_Matrix(fileLocation,fileName)

Loads just the S and T Matrices stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    Matrices = fload_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `SMatrix3` : A 6D matrix of the emission spectrum for 12->34 interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for 12->43 interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for 12->34 interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for 21->34 interaction.

If initial or final particles are identical then only one of the SMatrices or TMatrices will be returned for that state.

"""
function fload_Matrix(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        SMatrix3 = f["SMatrix3"];
        SMatrix4 = f["SMatrix4"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    # old format
    #(name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num) = Parameters

    if name1 == name2 && name3 == name4
        return (Parameters,SMatrix3,TMatrix1)
    end

    if name1 == name2 && name3 != name4
        return (Parameters,SMatrix3,SMatrix4,TMatrix1)
    end

    if name1 != name2 && name3 == name4
        return (Parameters,SMatrix3,TMatrix1,TMatrix2)
    end

    if name1 != name2 && name3 != name4
        return (Parameters,SMatrix3,SMatrix4,TMatrix1,TMatrix2)
    end

end

"""
    fload_Matrix_ISO(fileLocation,fileName)

Loads just the S and T Matrices stored in `fileName` stored at `fileLocation` first converting them to an isotropic form by summing over angles. (The dimensions of the matrices stay the same i.e. 6D->6D with three dimensions having a size of 1)

# Example
```julia-repl
    Matrices = fload_All_ISO(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Parameters` : A tuple of the parameters used in the evaluation.
- `SMatrix3` : A 6D matrix of the emission spectrum for 12->34 interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for 12->43 interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for 12->34 interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for 21->34 interaction.

If initial or final particles are identical then only one of the SMatrices or TMatrices will be returned for that state.

"""
function fload_Matrix_ISO(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        SMatrix3 = f["SMatrix3"];
        SMatrix4 = f["SMatrix4"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    if Parameters[1] == Parameters[2] && Parameters[3] == Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        return (Parameters,SMatrix3ISO,TMatrix1ISO)
    end

    if Parameters[1] == Parameters[2] && Parameters[3] != Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        SMatrix4ISO = sum(SMatrix4,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        return (Parameters,SMatrix3ISO,SMatrix4ISO,TMatrix1ISO)
    end

    if Parameters[1] != Parameters[2] && Parameters[3] == Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        TMatrix2ISO = sum(TMatrix2,dims=(2,4))
        return (Parameters,SMatrix3ISO,TMatrix1ISO,TMatrix2ISO)
    end

    if Parameters[1] != Parameters[2] && Parameters[3] != Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        SMatrix4ISO = sum(SMatrix4,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        TMatrix2ISO = sum(TMatrix2,dims=(2,4))
        return (SMatrix3ISO,SMatrix4ISO,TMatrix1ISO,TMatrix2ISO)
    end

end



"""
    DoesConserve(SMatrix3,SMatrix4,TMatrix1,TMatrix2,Parameters)

Function prints the ratio of the sum of the S and T matrices and their differences, for all interaction paths, as to check number and energy conservation for a particular interaction. Arguments are as outputted by the `fload_Matrix` function. 
"""
function DoesConserve2(Output::Tuple)

    # Output is tuple generated by fload_Matrix

    Parameters = Output[1] # always Output[1]

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Parameters

    if name1 == name2 && name3 == name4
        SMatrix3 = Output[2]
        SMatrix4 = zeros(size(SMatrix3))
        TMatrix1 = Output[3]
        TMatrix2 = zeros(size(TMatrix1))
    end

    if name1 == name2 && name3 != name4
        SMatrix3 = Output[2]
        SMatrix4 = Output[3]
        TMatrix1 = Output[4]
        TMatrix2 = zeros(size(TMatrix1))
    end

    if name1 != name2 && name3 == name4
        SMatrix3 = Output[2]
        SMatrix4 = zeros(size(SMatrix3))
        TMatrix1 = Output[3]
        TMatrix2 = Output[4]
    end

    if name1 != name2 && name3 != name4
        SMatrix3 = Output[2]
        SMatrix4 = Output[3]
        TMatrix1 = Output[4]
        TMatrix2 = Output[5]
    end

    #=mu1 = Parameters.mu1
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
    h1_num = Parameters.h1_num=#

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    p1_d_full = [p1_d; deltaVector([p1_r[end]; 2*p1_r[end]])];
    E1_Δ = deltaEVector(p1_r,mu1);
    E1_Δ_full = [E1_Δ; deltaEVector([p1_r[end], 2*p1_r[end]],mu1)];
    u1_r = bounds(u_low,u_up,u1_num,u1_grid);
    u1_d = deltaVector(u1_r);
    h1_r = bounds(h_low,h_up,h1_num,h1_grid) .* pi
    h1_d = deltaVector(h1_r);

    #=p2_low = Parameters.p2_low
    p2_up = Parameters.p2_up
    p2_grid = Parameters.p2_grid
    p2_num = Parameters.p2_num
    u2_grid = Parameters.u2_grid
    u2_num = Parameters.u2_num
    h2_grid = Parameters.h2_grid
    h2_num = Parameters.h2_num=#

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    p2_d_full = [p2_d; deltaVector([p2_r[end]; 2*p2_r[end]])];
    E2_Δ = deltaEVector(p2_r,mu2);
    E2_Δ_full = [E2_Δ; deltaEVector([p2_r[end], 2*p2_r[end]],mu2)];
    u2_r = bounds(u_low,u_up,u2_num,u2_grid);
    u2_d = deltaVector(u2_r);
    h2_r = bounds(h_low,h_up,h2_num,h2_grid) .* pi
    h2_d = deltaVector(h2_r);

    #=p3_low = Parameters.p3_low
    p3_up = Parameters.p3_up
    p3_grid = Parameters.p3_grid
    p3_num = Parameters.p3_num
    u3_grid = Parameters.u3_grid
    u3_num = Parameters.u3_num
    h3_grid = Parameters.h3_grid
    h3_num = Parameters.h3_num=#

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = deltaVector(p3_r);
    p3_d_full = [p3_d; deltaVector([p3_r[end]; 2*p3_r[end]])];
    E3_Δ = deltaEVector(p3_r,mu3);
    E3_Δ_full = [E3_Δ; deltaEVector([p3_r[end], 2*p3_r[end]],mu3)];
    u3_r = bounds(u_low,u_up,u3_num,u3_grid);
    u3_d = deltaVector(u3_r);
    h3_r = bounds(h_low,h_up,h3_num,h3_grid) .* pi
    h3_d = deltaVector(h3_r);

    #=p4_low = Parameters.p4_low
    p4_up = Parameters.p4_up
    p4_grid = Parameters.p4_grid
    p4_num = Parameters.p4_num
    u4_grid = Parameters.u4_grid
    u4_num = Parameters.u4_num
    h4_grid = Parameters.h4_grid
    h4_num = Parameters.h4_num=#

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_d = deltaVector(p4_r);
    p4_d_full = [p4_d; deltaVector([p4_r[end]; 2*p4_r[end]])];
    E4_Δ = deltaEVector(p4_r,mu4);
    E4_Δ_full = [E4_Δ; deltaEVector([p4_r[end], 2*p4_r[end]],mu4)];
    u4_r = bounds(u_low,u_up,u4_num,u4_grid);
    u4_d = deltaVector(u4_r);
    h4_r = bounds(h_low,h_up,h4_num,h4_grid) .* pi
    h4_d = deltaVector(h4_r);

    SsumN3 = 0
    TsumN1 = 0
    SsumE3 = 0
    TsumE1 = 0

    SsumN4 = 0
    TsumN2 = 0
    SsumE4 = 0
    TsumE2 = 0

    NGainMatrix3 = zeros(Float64,size(TMatrix1))
    NLossMatrix1 = zeros(Float64,size(TMatrix1))
    EGainMatrix3 = zeros(Float64,size(TMatrix1))
    ELossMatrix1 = zeros(Float64,size(TMatrix1))

    NGainMatrix4 = zeros(Float64,size(TMatrix1))
    NLossMatrix2 = zeros(Float64,size(TMatrix1))
    EGainMatrix4 = zeros(Float64,size(TMatrix1))
    ELossMatrix2 = zeros(Float64,size(TMatrix1))

    for p1 in axes(SMatrix3, 4), u1 in axes(SMatrix3,5), h1 in axes(SMatrix3,6), p2 in axes(SMatrix3,7), u2 in axes(SMatrix3,8), h2 in axes(SMatrix3,9)
        for p3 in axes(SMatrix3,1), u3 in axes(SMatrix3,2), h3 in axes(SMatrix3,3) 
        SsumN3 += SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*p3_d_full[p3]*u3_d[u3]*h3_d[h3]
        SsumE3 += SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_Δ_full[p3]*u3_d[u3]*h3_d[h3]
        NGainMatrix3[p1,u1,h1,p2,u2,h2] += SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*p3_d_full[p3]*u3_d[u3]*h3_d[h3]
        EGainMatrix3[p1,u1,h1,p2,u2,h2] += SMatrix3[p3,u3,h3,p1,u1,h1,p2,u2,h2]*E3_Δ_full[p3]*u3_d[u3]*h3_d[h3]
        end
    end

    for p1 in axes(SMatrix4, 4), u1 in axes(SMatrix4,5), h1 in axes(SMatrix4,6), p2 in axes(SMatrix4,7), u2 in axes(SMatrix4,8), h2 in axes(SMatrix4,9)
        for p4 in axes(SMatrix4,1), u4 in axes(SMatrix4,2), h4 in axes(SMatrix4,3) 
        SsumN4 += SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*p4_d_full[p4]*u4_d[u4]*h4_d[h4]
        SsumE4 += SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_Δ_full[p4]*u4_d[u4]*h4_d[h4]
        NGainMatrix4[p1,u1,h1,p2,u2,h2] += SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*p4_d_full[p4]*u4_d[u4]*h4_d[h4]
        EGainMatrix4[p1,u1,h1,p2,u2,h2] += SMatrix4[p4,u4,h4,p1,u1,h1,p2,u2,h2]*E4_Δ_full[p4]*u4_d[u4]*h4_d[h4]
        end
    end

    for p1 in axes(TMatrix1,1), u1 in axes(TMatrix1, 2), h1 in axes(TMatrix1,3), p2 in axes(TMatrix1,4), u2 in axes(TMatrix1,5), h2 in axes(TMatrix1,6)
        TsumN1 += TMatrix1[p1,u1,h1,p2,u2,h2]*p1_d_full[p1]*u1_d[u1]*h1_d[h1]
        TsumE1 += TMatrix1[p1,u1,h1,p2,u2,h2]*E1_Δ_full[p1]*u1_d[u1]*h1_d[h1]
        TsumN2 += TMatrix2[p2,u2,h2,p1,u1,h1]*p2_d_full[p2]*u2_d[u2]*h2_d[h2]
        TsumE2 += TMatrix2[p2,u2,h2,p1,u1,h1]*E2_Δ_full[p2]*u2_d[u2]*h2_d[h2]
        NLossMatrix1[p1,u1,h1,p2,u2,h2] += TMatrix1[p1,u1,h1,p2,u2,h2]*p1_d_full[p1]*u1_d[u1]*h1_d[h1]
        NLossMatrix2[p1,u1,h1,p2,u2,h2] += TMatrix2[p2,u2,h2,p1,u1,h1]*p2_d_full[p2]*u2_d[u2]*h2_d[h2]
        ELossMatrix1[p1,u1,h1,p2,u2,h2] += TMatrix1[p1,u1,h1,p2,u2,h2]*E1_Δ_full[p1]*u1_d[u1]*h1_d[h1]
        ELossMatrix2[p1,u1,h1,p2,u2,h2] += TMatrix2[p2,u2,h2,p1,u1,h1]*E2_Δ_full[p2]*u2_d[u2]*h2_d[h2]
    end

    NErrMatrix1 = (NGainMatrix3 .- NLossMatrix1) ./ NLossMatrix1
    EErrMatrix1 = (EGainMatrix3 .- ELossMatrix1) ./ ELossMatrix1
    meanNErr1 = sum(abs.(NErrMatrix1)) / length(NLossMatrix1)
    meanEErr1 = sum(abs.(EErrMatrix1)) / length(NLossMatrix1)
    stdN1 = sqrt(sum((NErrMatrix1 .- meanNErr1).^2)/length(NLossMatrix1))
    stdE1 = sqrt(sum((EErrMatrix1 .- meanEErr1).^2)/length(NLossMatrix1))

    NErrMatrix2 = (NGainMatrix4 .- NLossMatrix2) ./ NLossMatrix2
    EErrMatrix2 = (EGainMatrix4 .- ELossMatrix2) ./ ELossMatrix2
    meanNErr2 = sum(abs.(NErrMatrix2)) / length(NLossMatrix2)
    meanEErr2 = sum(abs.(EErrMatrix2)) / length(NLossMatrix2)
    stdN2 = sqrt(sum((NErrMatrix2 .- meanNErr2).^2)/length(NLossMatrix2))
    stdE2 = sqrt(sum((EErrMatrix2 .- meanEErr2).^2)/length(NLossMatrix2))

    println("sumSN3 = "*string(SsumN3))
    println("sumSN4 = "*string(SsumN4))
    println("sumTN1 = "*string(TsumN1))    
    println("sumTN2 = "*string(TsumN2)) 
    SsumN = SsumN3 + SsumN4
    println("sumSN = "*string(SsumN))
    TsumN = TsumN1 + TsumN2
    println("sumTN = "*string(TsumN))

    println("#")

    println("sumSE3 = "*string(SsumE3))
    println("sumSE4 = "*string(SsumE4))
    println("sumTE1 = "*string(TsumE1))  
    println("sumTE2 = "*string(TsumE2))
    SsumE = SsumE3 + SsumE4
    println("sumSE = "*string(SsumE))
    TsumE = TsumE1 + TsumE2
    println("sumTE = "*string(TsumE))

    println("#")

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string(SsumE-TsumE))
    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

    println("#")
    println("#")
    println("mean error in N = $meanNErr1")
    println("std of error in  N = $stdN1")
    println("mean error in E = $meanEErr1")
    println("std of error in E = $stdE1")
    println("#")
    println("#")
    println("mean error in N = $meanNErr2")
    println("std of error in  N = $stdN2")
    println("mean error in E = $meanEErr2")
    println("std of error in E = $stdE2")

    return NGainMatrix3, NLossMatrix1, NErrMatrix1,NGainMatrix4, NLossMatrix2, NErrMatrix2

end


"""
    DoesConserve(SMatrix3,SMatrix4,TMatrix1,TMatrix2,Parameters)

Function prints the ratio of the sum of the S and T matrices and their differences, for all interaction paths, as to check number and energy conservation for a particular interaction. Arguments are as outputted by the `fload_All` function. 
"""
function DoesConserve(SMatrix3,SMatrix4,TMatrix1,TMatrix2,Parameters)

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num) = Parameters

    mu1 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
    mu2 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
    mu3 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
    mu4 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

    p1_r = bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = deltaVector(p1_r);
    p1_d_full = [p1_d; deltaVector([p1_r[end]; 2*p1_r[end]])];
    E1_Δ = deltaEVector(p1_r,mu1);
    E1_Δ_full = [E1_Δ; deltaEVector([p1_r[end], 2*p1_r[end]],mu1)];
    u1_r = bounds(u_low,u_up,u1_num,u1_grid);
    u1_d = deltaVector(u1_r);

    p2_r = bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = deltaVector(p2_r);
    p2_d_full = [p2_d; deltaVector([p2_r[end]; 2*p2_r[end]])];
    E2_Δ = deltaEVector(p2_r,mu2);
    E2_Δ_full = [E2_Δ; deltaEVector([p2_r[end], 2*p2_r[end]],mu2)];
    u2_r = bounds(u_low,u_up,u2_num,u2_grid);
    u2_d = deltaVector(u2_r);

    p3_r = bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = deltaVector(p3_r);
    p3_d_full = [p3_d; deltaVector([p3_r[end]; 2*p3_r[end]])];
    E3_Δ = deltaEVector(p3_r,mu3);
    E3_Δ_full = [E3_Δ; deltaEVector([p3_r[end], 2*p3_r[end]],mu3)];
    u3_r = bounds(u_low,u_up,u3_num,u3_grid);
    u3_d = deltaVector(u3_r);

    p4_r = bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_d = deltaVector(p4_r);
    p4_d_full = [p4_d; deltaVector([p4_r[end]; 2*p4_r[end]])];
    E4_Δ = deltaEVector(p4_r,mu4);
    E4_Δ_full = [E4_Δ; deltaEVector([p4_r[end], 2*p4_r[end]],mu4)];
    u4_r = bounds(u_low,u_up,u4_num,u4_grid);
    u4_d = deltaVector(u4_r);

    SsumN3 = 0
    TsumN1 = 0
    SsumE3 = 0
    TsumE1 = 0

    SsumN4 = 0
    TsumN2 = 0
    SsumE4 = 0
    TsumE2 = 0

    for k in axes(SMatrix3,3), l in axes(SMatrix3, 4), m in axes(SMatrix3,5), n in axes(SMatrix3,6)
        for i in axes(SMatrix3,1), j in axes(SMatrix3,2) 
        SsumN3 += SMatrix3[i,j,k,l,m,n]*p3_d_full[i]*u3_d[j]
        SsumE3 += SMatrix3[i,j,k,l,m,n]*E3_Δ_full[i]*u3_d[j]
        end
    end

    for k in axes(SMatrix4,3), l in axes(SMatrix4, 4), m in axes(SMatrix4,5), n in axes(SMatrix4,6)
        for i in axes(SMatrix4,1), j in axes(SMatrix4,2) 
        SsumN4 += SMatrix4[i,j,k,l,m,n]*p4_d_full[i]*u4_d[j]
        SsumE4 += SMatrix4[i,j,k,l,m,n]*E4_Δ_full[i]*u4_d[j]
        end
    end

    for k in axes(TMatrix1,1), l in axes(TMatrix1, 2), m in axes(TMatrix1,3), n in axes(TMatrix1,4)
        TsumN1 += TMatrix1[k,l,m,n]*p1_d_full[k]*u1_d[l]
        TsumE1 += TMatrix1[k,l,m,n]*E1_Δ_full[k]*u1_d[l]
        TsumN2 += TMatrix2[m,n,k,l]*p2_d_full[m]*u2_d[n]
        TsumE2 += TMatrix2[m,n,k,l]*E2_Δ_full[m]*u2_d[n]
    end

    println("sumSN3 = "*string(SsumN3))
    println("sumSN4 = "*string(SsumN4))
    println("sumTN1 = "*string(TsumN1))    
    println("sumTN2 = "*string(TsumN2)) 
    SsumN = SsumN3 + SsumN4
    println("sumSN = "*string(SsumN))
    TsumN = TsumN1 + TsumN2
    println("sumTN = "*string(TsumN))

    println("#")

    println("sumSE3 = "*string(SsumE3))
    println("sumSE4 = "*string(SsumE4))
    println("sumTE1 = "*string(TsumE1))  
    println("sumTE2 = "*string(TsumE2))
    SsumE = SsumE3 + SsumE4
    println("sumSE = "*string(SsumE))
    TsumE = TsumE1 + TsumE2
    println("sumTE = "*string(TsumE))

    println("#")

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string((SsumE-TsumE)))
    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

end


