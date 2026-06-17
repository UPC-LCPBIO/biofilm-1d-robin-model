# Compatibility entry point for the paper reproduction scripts.
# Including this file exposes the reduced-model solver and helper modules in the
# caller namespace, matching the scripts used to generate the manuscript data.

include(joinpath(@__DIR__, "solver.jl"))
