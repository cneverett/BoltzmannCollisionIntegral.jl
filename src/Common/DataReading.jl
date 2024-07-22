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