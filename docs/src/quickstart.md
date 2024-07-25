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
- Define the names of the 4 particles involved in the interaction (12->34) as the strings `name1` `name2` `name3` `name4`
    - These should of the form of three letters, which abbreviate the particles full name (see [Particles](@ref) for list of currently implemented particles).
    - They should be ordered to match a currently [Implemented Interactions](@ref)
- Define the momentum space discretisation. This includes the upper and lower bounds of momentum for particle species 1,2 and 3 (e.g. `p1l` and `p1u` for species 1) and the number of divisions (bins) for each particles momentum space (e.g. `nump1`).
- Define the number of divisions for the angular momentum space (i.e. cos(theta) space) for the particle species 1,2 and 3 (e.g. `numt1`). 
- Define the number of Monte-Carlo samples to perform (as a rule of thumb, on a modern CPU it takes approximately 200ns per sample).
    - `numTiter` for the number of random sets of $\{\vec{p}_1,\vec{p}_2\}$ to sample. 
    - `numSiter` for the number of random $\{\vec{p}_3\}$ states to sample per $\{\vec{p}_1,\vec{p}_2\}$.
- If multithreading then define `nThreads` that will be used. This generates `nThreads` workers that perform evaluation in parallel, utilising `locks` to prevent data races. (see [Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/) for how to set up multi-threading in Julia)
- Define the `fileLocation` where the output file ([JLD2](https://github.com/JuliaIO/JLD2.jl)) named `fileName` is to be written.
- Evaluate the emission and absorption spectrum using the [`SpectraEvaluateSerial`](@ref) function for serial and [SpectraEvaluateMultiThread](@ref) for multithread. Once run, these functions will save the results to the output file.

## Output
- The data stored in an output file can be loaded back into the workspace as a tuple using the [fload_All](@ref) function.
- Conservation of particle number and energy can be checked using the [DoesConserve](@ref) function.

