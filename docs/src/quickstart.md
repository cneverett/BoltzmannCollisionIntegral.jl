# Getting Started

First install the package by 
... TO BE UPDATED ONCE PACKAGE IS PUBLICLY AVAILABLE

!!! note
    An example script `Run_Integration.jl` for setting up and running the evaluation of the discrete collision integral can be found under `src/Common/` of the package. It is recommended that you copy this script to your working directory, edit the relevant fields and then running, either using the command 
    ```julia-repl
    include("Run_Integration.jl)
    ```
    in a julia-repl session, or by running the script line by line in your favourite code editor.

The example script `Run_Integration.jl` operates as follows:
- Define the names of the 4 particles involved in the interaction (``12\to34``) as the strings `name1` `name2` `name3` `name4`
    - These should of the form of three letters, which abbreviate the particles full name (see [Particles](@ref) for list of currently implemented particles).
    - They should be ordered to match a currently [Implemented Interactions](@ref)
- Define the momentum space discretisation. This includes the upper and lower bounds of momentum (log10) for particle species 1,2 and 3 (e.g. `p1l` and `p1u` for species 1) and the number of divisions (bins) for each particles momentum space (e.g. `nump1`).
- Define the number of divisions for the angular momentum space (i.e. cos(theta) space) for the particle species 1,2 and 3 (e.g. `numt1`). 
- Define the number of Monte-Carlo samples to perform (as a rule of thumb, on a modern CPU it takes approximately 200ns per sample).
    - `numTiter` for the number of random sets of $\{\vec{p}_1,\vec{p}_2\}$ to sample. 
    - `numSiter` for the number of random $\{\vec{p}_3\}$ states to sample per $\{\vec{p}_1,\vec{p}_2\}$.
- If multithreading then define `nThreads` that will be used. This generates `nThreads` workers that perform evaluation in parallel, utilising `locks` to prevent data races. (see [Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to set up multi-threading in Julia)
- Define the `fileLocation` where the output file ([JLD2](https://github.com/JuliaIO/JLD2.jl)) named `fileName` is to be written.
- Evaluate the emission and absorption spectrum using the [`SpectraEvaluateSerial`](@ref) function for serial and [SpectraEvaluateMultiThread](@ref) for multithread. Once run, these functions will save the results to the output file.

## Output
The data stored in an output file can be loaded back into the workspace as a tuple using the [fload_All](@ref) function.
        
```julia-repl
    (Run_Parameters,SAtot3,SAtot4,TAtot,SAtal3,SAtal4,TAtal,SMatrix3,SMatrix4,TMatrix1,TMatrix2,p3Max,p4Max,t3MinMax,t4MinMax,SConv3,SConv4,TConv) = fload_All(fileLocation,fileName)
```

This returns a `Tuple` containing various arrays:
- `Run_Parameters` : A tuple of the parameters used in the evaluation.
- `Stot3` : A 6D matrix totalling all the emission spectrum values sampled for ``12\to34`` interaction.
- `Stot4` : A 6D matrix totalling all the emission spectrum values sampled for ``12\to43`` interaction.
- `Ttot` : A 4D matrix totalling all the absorption spectrum values sampled.
- `Stal3` : A 5D matrix of tallies of the number of emission spectrum values sampled for 12->34 interaction.
- `Stal4` : A 5D matrix of tallies of the number of emission spectrum values sampled for 12->43 interaction.
- `Ttal` : A 4D matrix of tallies of the number of absorption spectrum values sampled.
- `SMatrix3` : A 6D matrix of the emission spectrum for ``12\to34`` interaction.
- `SMatrix4` : A 6D matrix of the emission spectrum for ``12\to43`` interaction.
- `TMatrix1` : A 4D matrix of the absorption spectrum for ``12\to34`` interaction.
- `TMatrix2` : A 4D matrix of the absorption spectrum for ``21\to34`` interaction i.e. by permutation of TMatrix1 and correct application of phase space factors if species 1 != species 2.
- `p3Max` : The maximum value of the momentum space variable p3 sampled for each bin. (Useful for correcting numerical diffusion)
- `t3MinMax` : The minimum and maximum values of the momentum space variable t3 sampled for each bin. (Useful for correcting numerical diffusion)
- `p4Max` : The maximum value of the momentum space variable p4 sampled for each bin. (Useful for correcting numerical diffusion)
- `t4MinMax` : The minimum and maximum values of the momentum space variable t4 sampled for each bin. (Useful for correcting numerical diffusion)
- `SConv3` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Run_Parameters` for ``12\to34`` interaction.
- `SConv4` : A 6D matrix of the convergence of the emission spectrum compared to the previous run with given `Run_Parameters` for ``12\to43`` interaction.
- `TConv` : A 4D matrix of the convergence of the absorption spectrum compared to the previous run with given `Run_Parameters`.

Conservation of particle number and energy can be checked using the [DoesConserve](@ref) function.
- The key statistic is `ratioN` and `ratioE` which dictate the ratio of particle number and energy before and after the interaction and should be close to 1.

