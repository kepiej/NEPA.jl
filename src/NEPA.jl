__precompile__(true)
module NEPA
# Nonparametric Efficiency and Productivity Analysis (NEPA)

using MathProgBase
using Clp

# Export functions and types
export RS, CRS, VRS, NIRS, NDRS
export AbstractDEA, AbstractDataEnvelopment, getdata
export DEAData, getnrdmu, getiodim, getindexes, setindexes!, size, linearindexing, getindex
export DEAResult, geteff, getpeers, getsx, getsy, show
export Convex, FreeDisposal, DDF
export DEA, DEA_CRS, DEA_VRS, DEA_NIRS, DEA_NDRS
export FDH, FDH_CRS, FDH_VRS, FDH_NIRS, FDH_NDRS
export SBM
export WACM
export Luenberger,TEI,TC
export HMB,TEI,TC,SEC
export orderm
export Hyperbolic

# Source files
include("RTS.jl")
include("AbstractDEA.jl")
include("DEAData.jl")
include("DEAResult.jl")
include("DDF.jl")
include("DEA.jl")
include("FDH.jl")
include("SBM.jl")
include("WACM.jl")
include("Luenberger.jl")
include("HMB.jl")
include("OrderM.jl")
include("Hyperbolic.jl")

end
