module ADI

using ConcreteStructs
using LinearAlgebra
using ProgressLogging
using Reexport
using Statistics
using StructArrays

@reexport using HCIToolbox

export reconstruct,
       subtract,
       process,
       Classic,
       LOCI,
       PCA,
       GreeDS,
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
include("loci.jl")
include("pca.jl")
include("nmf.jl")
include("greeds.jl")

# metrics
include("metrics/Metrics.jl")
@reexport using .Metrics

end
