"""
    fload_All_Sync(fileLocation,fileName)

Loads all the data stored in `fileName` stored at `fileLocation`.

# Example
```julia-repl
    (Run_Parameters,SAtot,SAtal,SMatrix,#=pMax,tMinMax,=#SConv) = fload_All_Sync(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `Run_Parameters` : A tuple of the parameters used in the evaluation.
- `Stot` : A 4D matrix totalling all the synchrotron emission spectrum values
- `Stal` : A 4D matrix of tallies of the number of synchrotron emission spectrum values sampled
- `SMatrix` : A 4D matrix of the synchrotron emission spectrum.
- `pMax` : The maximum value of the momentum space variable p1 (photon mommentum) sampled for each bin. (Useful for correcting numerical diffusion)
- `tMinMax` : The minimum and maximum values of the momentum space variable t1 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv` : A 4D matrix of the convergence of the synchrotron emission spectrum compared to the previous run with given `Run_Parameters`.
"""
function fload_All_Sync(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Run_Parameters = f["Parameters"]

        SAtot = f["STotal"];
        SAtal = f["STally"];
        SMatrix = f["SMatrix"];
        SConv = f["SConverge"];
        #pMax = f["pMax"];
        #tMinMax = f["tMinMax"];

        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Run_Parameters,SAtot,SAtal,SMatrix,#=pMax,tMinMax,=#SConv);

end