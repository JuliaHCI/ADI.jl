using ADI
using Documenter
using HCIToolbox
using Literate: markdown

# weave in examples using literate
examplesdir(args...) = joinpath(@__DIR__, "..", "examples", args...)
outdir = joinpath(@__DIR__, "src", "examples")
markdown(examplesdir("betapictoris.jl"), outdir)
markdown(examplesdir("geometries.jl"), outdir)

# now do Documenter commands

setup = quote
    using ADI
    using Random
    Random.seed!(8799)
end

DocMeta.setdocmeta!(ADI, :DocTestSetup, setup; recursive = true)

doctest(ADI)

makedocs(;
    modules = [ADI, HCIToolbox],
    authors = "Miles Lucas <mdlucas@hawaii.edu>",
    repo = "https://github.com/juliahci/ADI.jl/blob/{commit}{path}#L{line}",
    sitename = "ADI.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliahci.github.io/ADI.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Introduction to HCI" => "introduction.md",
        "Getting Started with ADI.jl" => "gettingstarted.md",
        "Algorithms" => [
            "algorithms/classic.md",
            "algorithms/pca.md",
            "algorithms/nmf.md",
            "algorithms/greeds.md",
            "algorithms/api.md"
        ],
        "Framewise" => "framewise.md",
        "SDI" => "sdi.md",
        "Metrics" => "metrics.md",
        "Examples" => [
            "examples/betapictoris.md",
            "examples/geometries.md"
        ],
        "Benchmarks" => "benchmarks.md",
        "Index" => "api.md"
    ],
    doctest = false
)

deploydocs(;
    repo = "github.com/JuliaHCI/ADI.jl",
    push_preview = true
)
