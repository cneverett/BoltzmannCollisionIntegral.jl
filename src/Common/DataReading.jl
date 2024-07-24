"""
    fload_All(fileLocation,fileName)

Loads all the data stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    (Run_Parameters, Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Run_Parameters` : A tuple of the parameters used in the evaluation.
- `Stot` : A 6D matrix totalling all the emission spectrum values sampled.
- `Ttot` : A 4D matrix totalling all the absorption spectrum values sampled.
- `Stal` : A 5D matrix of tallies of the number of emission spectrum values sampled.
- `Ttal` : A 4D matrix of tallies of the number of absorption spectrum values sampled.
- `SMatrix` : A 6D matrix of the emission spectrum.
- `TMatrix` : A 4D matrix of the absorption spectrum.
- `p3Max` : The maximum value of the momentum space variable p3 sampled for each bin. (Useful for correcting numerical diffusion)
- `t3MinMax` : The minimum and maximum values of the momentum space variable t3 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv` : A 6D matrix of the convergence of the emission spectrum compaired to the previous run with given `Run_Parameters`.
- `TConv` : A 4D matrix of the convergence of the absorption spectrum compaired to the previous run with given `Run_Parameters`.

"""
function fload_All(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Run_Parameters = f["Parameters"]
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
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Run_Parameters, SAtot, TAtot, AStal, ATtal, SMatrix, TMatrix, p3Max, t3MinMax, SConv, TConv);
    #run (Run_Parameters, Stot,Ttot,Stal,Ttal,SMatrix,TMatrix,p3Max,t3MinMax,SConv,TConv) = fload_All(fileLocation,fileName); in REPL

end

"""
    DoesConserve(SMatrix,TMatrix,Parameters)

Function prints the ratio of the sum of the S and T matricies and their differences as to check number and energy conservation for a particular interaction. Arguments are as outputted by the `fload_All` function. 
"""
function DoesConserve(SMatrix,TMatrix,Run_Parameters)

    (name1,name2,name3,name4,p3l,p3u,nump3,p1l,p1u,nump1,p2l,p2u,nump2,numt3,numt1,numt2) = Run_Parameters

    mu1 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name1))
    mu2 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name2))
    mu3 = getfield(BoltzmannCollisionIntegral,Symbol("mu"*name3))

    pr3 = prange(p3l,p3u,nump3);
    dp3 = deltaVector(pr3);
    dp3full = [dp3; deltaVector([pr3[end]; 2*pr3[end]])];
    ΔE3 = deltaEVector(pr3,mu3);
    ΔE3full = [ΔE3; deltaEVector([pr3[end], 2*pr3[end]],mu3)];
    tr3 = trange(numt3);
    dμ3 = deltaVector(tr3);

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

    SsumN = 0f0
    TsumN = 0f0
    SsumE = 0f0
    TsumE = 0f0

    for k in axes(SMatrix,3), l in axes(SMatrix, 4), m in axes(SMatrix,5), n in axes(SMatrix,6)
        for i in axes(SMatrix,1), j in axes(SMatrix,2) 
        SsumN += SMatrix[i,j,k,l,m,n]*dp3full[i]*dμ3[j]
        SsumE += SMatrix[i,j,k,l,m,n]*ΔE3full[i]*dμ3[j]
        end
    end

    for k in axes(TMatrix,1), l in axes(TMatrix, 2), m in axes(TMatrix,3), n in axes(TMatrix,4)
        TsumN += TMatrix[k,l,m,n]*dp1full[k]*dμ1[l]
        TsumE += TMatrix[k,l,m,n]*ΔE1full[k]*dμ1[l]
        if name1!=name2 # total absorption if particles are not the same
        #    TsumN += TMatrix[k,l,m,n]*dp2full[m]*dμ2[n]
        #    TsumE += TMatrix[k,l,m,n]*ΔE2full[m]*dμ2[n]
        end
    end

    println("sumSN = "*string(SsumN))
    println("sumTN = "*string(TsumN))    

    println("sumSE = "*string(SsumE))
    println("sumTE = "*string(TsumE))  

    println("errN = "*string(SsumN-TsumN))
    println("errE = "*string((SsumE-TsumE)))

    println("ratioN = "*string(SsumN/TsumN))
    println("ratioE = "*string(SsumE/TsumE))

end