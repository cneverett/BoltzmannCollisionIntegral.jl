# Inside make.jl
push!(LOAD_PATH,"../src/")
using BinaryInteractionSpectra
using Documenter
makedocs(
         sitename = "BinaryInteractionSpectra.jl",
         modules  = [BinaryInteractionSpectra],
         pages=[
                "Overview" => "index.md",
                "Quick Start" => "quickstart.md",
                "Cross Sections" => "crosssections.md"
                "Internal Functions" => "internalfunctions.md"
               ])
deploydocs(;
    repo="github.com/cneverett/BinaryInteractionSpectra.jl",
)