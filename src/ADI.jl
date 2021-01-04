module ADI

using LinearAlgebra
using ProgressMeter
using Reexport
using Statistics
using StructArrays

@reexport using HCIToolbox

export reconstruct,
       subtract,
       decompose,
       Classic,
       PCA,
       GreeDS,
       TPCA,
       NMF,
       Framewise,
       SingleSDI,
       DoubleSDI,
       SliceSDI

# core API and default implementations
include("common.jl")
include("design.jl")

# further techniques
include("framewise.jl")
include("sdi.jl")

# algorithms
include("classic.jl")
include("pca.jl")
include("nmf.jl")
include("greeds.jl")

# metrics
include("metrics/Metrics.jl")
@reexport using .Metrics

end
