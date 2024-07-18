# Inside make.jl
push!(LOAD_PATH,"../src/")
using BinaryInteractionSpectra
using Documenter
makedocs(
         sitename = "BinaryInteractionSpectra.jl",
         modules  = [BinaryInteractionSpectra],
         pages=[
                "Overview" => "index.md",
                "Getting Started" => "quickstart.md",
                "Cross Sections" => "crosssections.md",
                "Internal Functions" => "internalfunctions.md"
               ],
        #format = Documenter.HTML(
        #    mathengine = MathJax3(Dict(
        #        :loader => Dict("load" => ["[tex]/physics"]),
        #        :tex => Dict(
        #            "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        #            "tags" => "ams",
        #            "packages" => ["base", "ams", "autoload", "physics"],
        #        ),
        #    ))
        #),
        )
deploydocs(;
    repo="github.com/cneverett/BinaryInteractionSpectra.jl",
)