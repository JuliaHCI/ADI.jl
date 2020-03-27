using ADI
using Documenter

setup = quote
    using ADI
    using Random
    Random.seed!(8799)
end

DocMeta.setdocmeta!(ADI, :DocTestSetup, setup; recursive = true)


makedocs(;
    modules = [ADI],
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
        "API/Reference" => "api.md"
    ],
)

deploydocs(;
    repo = "github.com/JuliaHCI/ADI.jl",
    push_preview = true
)
