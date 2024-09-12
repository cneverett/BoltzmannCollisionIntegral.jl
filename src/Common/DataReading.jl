"""
    fload_All(fileLocation,fileName)

Loads all the data stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    (Run_Parameters,SAtot3,SAtot4,TAtot,SAtal3,SAtal4,TAtal,SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3Max,p4Max,t3MinMax,t4MinMax,SConv3,SConv4,TConv) = fload_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Run_Parameters` : A tuple of the parameters used in the evaluation.
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
- `t3MinMax` : The minimum and maximum values of the momentum space variable t3 sampled for each bin. (Useful for correcting numerical diffusion)
- `p4Max` : The maximum value of the momentum space variable p4 sampled for each bin. (Useful for correcting numerical diffusion)
- `t4MinMax` : The minimum and maximum values of the momentum space variable t4 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv3` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Run_Parameters` for 12->34 interaction.
- `SConv4` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Run_Parameters` for 12->43 interaction.
- `TConv` : A 4D matrix of the convergence of the absorption spectrum compared to the previous run with given `Run_Parameters`.

"""
function fload_All(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Run_Parameters = f["Parameters"]

        SAtot3 = f["STotal3"];
        SAtal3 = f["STally3"];
        SMatrix3 = f["SMatrix3"];
        SConv3 = f["SConverge3"];
        p3Max = f["p3Max"];
        t3MinMax = f["t3MinMax"];

        SAtot4 = f["STotal4"];
        SAtal4 = f["STally4"];
        SMatrix4 = f["SMatrix4"];
        SConv4 = f["SConverge4"];
        p4Max = f["p4Max"];
        t4MinMax = f["t4MinMax"];

        TAtot = f["TTotal"];
        TAtal = f["TTally"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];
        TConv = f["TConverge"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Run_Parameters,SAtot3,SAtot4,TAtot,SAtal3,SAtal4,TAtal,SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3Max,p4Max,t3MinMax,t4MinMax,SConv3,SConv4,TConv);

end


"""
    fload_Matrix(fileLocation,fileName)

Loads just the S and T Matricies stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    Matricies = fload_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `SMatrix3` : A 6D matrix of the emission spectrum for 12->34 interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for 12->43 interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for 12->34 interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for 21->34 interaction.

If initial or final particles are identical then only one of the SMatricies or TMatricies will be returned for that state.

"""
function fload_Matrix(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Run_Parameters = f["Parameters"]
        SMatrix3 = f["SMatrix3"];
        SMatrix4 = f["SMatrix4"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    if Run_Parameters[1] == Run_Parameters[2] && Run_Parameters[3] == Run_Parameters[4]
        return (SMatrix3,TMatrix1)
    end

    if Run_Parameters[1] == Run_Parameters[2] && Run_Parameters[3] != Run_Parameters[4]
        return (SMatrix3,SMatrix4,TMatrix1)
    end

    if Run_Parameters[1] != Run_Parameters[2] && Run_Parameters[3] == Run_Parameters[4]
        return (SMatrix3,TMatrix1,TMatrix2)
    end

    if Run_Parameters[1] != Run_Parameters[2] && Run_Parameters[3] != Run_Parameters[4]
        return (SMatrix3,SMatrix4,TMatrix1,TMatrix2)
    end

end

"""
    fload_Matrix_ISO(fileLocation,fileName)

Loads just the S and T Matricies stored in `fileName` stored at `fileLocation` first converting them to an isotropic form by summing over angles. (The dimensions of the matricies stay the same i.e. 6D->6D with three dimensions having a size of 1)

# Example
```julia-repl
    Matricies = fload_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `SMatrix3` : A 6D matrix of the emission spectrum for 12->34 interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for 12->43 interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for 12->34 interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for 21->34 interaction.

If initial or final particles are identical then only one of the SMatricies or TMatricies will be returned for that state.

"""
function fload_Matrix_ISO(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Run_Parameters = f["Parameters"]
        SMatrix3 = f["SMatrix3"];
        SMatrix4 = f["SMatrix4"];
        TMatrix1 = f["TMatrix1"];
        TMatrix2 = f["TMatrix2"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    if Run_Parameters[1] == Run_Parameters[2] && Run_Parameters[3] == Run_Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        return (SMatrix3ISO,TMatrix1ISO)
    end

    if Run_Parameters[1] == Run_Parameters[2] && Run_Parameters[3] != Run_Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        SMatrix4ISO = sum(SMatrix4,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        return (SMatrix3ISO,SMatrix4ISO,TMatrix1ISO)
    end

    if Run_Parameters[1] != Run_Parameters[2] && Run_Parameters[3] == Run_Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        TMatrix2ISO = sum(TMatrix2,dims=(2,4))
        return (SMatrix3ISO,TMatrix1ISO,TMatrix2ISO)
    end

    if Run_Parameters[1] != Run_Parameters[2] && Run_Parameters[3] != Run_Parameters[4]
        SMatrix3ISO = sum(SMatrix3,dims=(2,4,6))
        SMatrix4ISO = sum(SMatrix4,dims=(2,4,6))
        TMatrix1ISO = sum(TMatrix1,dims=(2,4))
        TMatrix2ISO = sum(TMatrix2,dims=(2,4))
        return (SMatrix3ISO,SMatrix4ISO,TMatrix1ISO,TMatrix2ISO)
    end

end


"""
    DoesConserve(SMatrix,TMatrix,Parameters)

Function prints the ratio of the sum of the S and T matricies and their differences, for all interaction paths, as to check number and energy conservation for a particular interaction. Arguments are as outputted by the `fload_All` function. 
"""
function DoesConserve(SMatrix3,SMatrix4,TMatrix1,TMatrix2,Run_Parameters)

    (name1,name2,name3,name4,p3l,p3u,nump3,p4l,p4u,nump4,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt4,numt1,numt2) = Run_Parameters

    mu1 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
    mu2 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
    mu3 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))
    mu4 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name4))

    pr3 = prange(p3l,p3u,nump3);
    dp3 = deltaVector(pr3);
    dp3full = [dp3; deltaVector([pr3[end]; 2*pr3[end]])];
    ΔE3 = deltaEVector(pr3,mu3);
    ΔE3full = [ΔE3; deltaEVector([pr3[end], 2*pr3[end]],mu3)];
    tr3 = trange(numt3);
    dμ3 = deltaVector(tr3);

    pr4 = prange(p4l,p4u,nump4);
    dp4 = deltaVector(pr4);
    dp4full = [dp4; deltaVector([pr4[end]; 2*pr4[end]])];
    ΔE4 = deltaEVector(pr4,mu4);
    ΔE4full = [ΔE4; deltaEVector([pr4[end], 2*pr4[end]],mu4)];
    tr4 = trange(numt4);
    dμ4 = deltaVector(tr4);

    pr1 = prange(p1l,p1u,nump1);
    dp1 = deltaVector(pr1);
    dp1full = [dp1; deltaVector([pr1[end]; 2*pr1[end]])];
    ΔE1 = deltaEVector(pr1,mu1);
    ΔE1full = [ΔE1; deltaEVector([pr1[end], 2*pr1[end]],mu1)];
    tr1 = trange(numt1);
    dμ1 = deltaVector(tr1);

    pr2 = prange(p2l,p2u,nump2);
    dp2 = deltaVector(pr2);
    dp2full = [dp2; deltaVector([pr2[end]; 2*pr2[end]])];
    ΔE2 = deltaEVector(pr2,mu2);
    ΔE2full = [ΔE2; deltaEVector([pr2[end], 2*pr2[end]],mu2)];
    tr2 = trange(numt2);
    dμ2 = deltaVector(tr2);

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
        SsumN3 += SMatrix3[i,j,k,l,m,n]*dp3full[i]*dμ3[j]
        SsumE3 += SMatrix3[i,j,k,l,m,n]*ΔE3full[i]*dμ3[j]
        end
    end

    for k in axes(SMatrix4,3), l in axes(SMatrix4, 4), m in axes(SMatrix4,5), n in axes(SMatrix4,6)
        for i in axes(SMatrix4,1), j in axes(SMatrix4,2) 
        SsumN4 += SMatrix4[i,j,k,l,m,n]*dp4full[i]*dμ4[j]
        SsumE4 += SMatrix4[i,j,k,l,m,n]*ΔE4full[i]*dμ4[j]
        end
    end

    for k in axes(TMatrix1,1), l in axes(TMatrix1, 2), m in axes(TMatrix1,3), n in axes(TMatrix1,4)
        TsumN1 += TMatrix1[k,l,m,n]*dp1full[k]*dμ1[l]
        TsumE1 += TMatrix1[k,l,m,n]*ΔE1full[k]*dμ1[l]
        TsumN2 += TMatrix2[m,n,k,l]*dp2full[m]*dμ2[n]
        TsumE2 += TMatrix2[m,n,k,l]*ΔE2full[m]*dμ2[n]
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


