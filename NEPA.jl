__precompile__()
module NEPA
# Nonparametric Efficiency and Productivity Analysis (NEPA)

using MathProgBase

# Export functions and types
export RS, CRS, VRS, NIRS, NDRS
export DEA, DEA_CRS, DEA_VRS, DEA_NIRS, DEA_NDRS
export FDH, FDH_CRS, FDH_VRS, FDH_NIRS, FDH_NDRS
export WACM

# Source files
include("RTS.jl")
include("DEA.jl")
include("FDH.jl")
include("WACM.jl")

end