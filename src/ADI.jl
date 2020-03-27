module ADI

using HCIToolbox
import MultivariateStats

const mvs = MultivariateStats

export pca

# The core decomposition routines
include("pca.jl")

end
