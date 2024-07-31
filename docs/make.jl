# Inside make.jl
push!(LOAD_PATH,"../src/")
using BoltzmannCollisionIntegral
using Documenter

makedocs(
        sitename = "BoltzmannCollisionIntegral.jl",
        modules  = [BoltzmannCollisionIntegral],
        authors = "Christopher Everett",
        pages=[
                "Overview" => "index.md",
                "Getting Started" => "quickstart.md",
                "Cross Sections" => "crosssections.md",
                "Internal Functions" => "internalfunctions.md"
               ],
        checkdocs = :export
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

deploydocs(
    target = "build",
    repo="github.com/cneverett/BoltzmannCollisionIntegral.jl.git",
    branch = "gh-pages",
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"]
)