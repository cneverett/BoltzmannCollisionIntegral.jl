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
        Parameters = f["Parameters"]

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

    return (Parameters,SAtot,SAtal,SMatrix,#=pMax,tMinMax,=#SConv);

end

"""
    fload_Matrix_Sync(fileLocation,fileName)

Loads just the S and T Matrices stored in `fileName` stored at `fileLocation`. 

# Example
```julia-repl
    Matrices = fload_Matrix_Sync(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `SMatrix` : A 4D matrix of the emission spectrum for Synchrotron.

"""
function fload_Matrix_Sync(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        GainMatrix3 = f["GainMatrix3"];
        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    return (Parameters,GainMatrix3)

end

"""
    fload_Matrix_SyncISO(fileLocation,fileName)

Loads just the S and T Matrices stored in `fileName` stored at `fileLocation` first converting them to an isotropic form by summing over angles. (The dimensions of the matrices stay the same i.e. 6D->6D with three dimensions having a size of 1)

# Example
```julia-repl
    Matrices = fload_Matrix_SyncISO(fileLocation,fileName);
```
Returns a tuple of the data stored in the file. The fields are as follows:
- `SMatrix` : A 4D matrix of the emission spectrum for Synchrotron.

"""
function fload_Matrix_SyncISO(fileLocation::String,fileName::String)
        
    filePath = fileLocation*"\\"*fileName
    fileExist = isfile(filePath)

    if fileExist
        f = jldopen(filePath,"r+");
        Parameters = f["Parameters"]
        SMatrix = f["SMatrix"];
        close(f)  
    else
        error("no file with name $fileName found at location $fileLocation")
    end

    SMatrixISO = sum(SMatrix,dims=(2,4))
    return (Parameters,SMatrixISO)

end