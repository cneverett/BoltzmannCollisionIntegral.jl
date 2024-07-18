# Inside make.jl
push!(LOAD_PATH,"../src/")
using BinaryInteractionSpectra
using Documenter
makedocs(
         sitename = "BinaryInteractionSpectra.jl",
         modules  = [BinaryInteractionSpectra],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/cneverett/BinaryInteractionSpectra.jl",
)