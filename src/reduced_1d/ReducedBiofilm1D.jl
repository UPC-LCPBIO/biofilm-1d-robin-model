module ReducedBiofilm1D

include("solver.jl")
include("hydrodynamics.jl")
include("campaigns.jl")

export run_case
export threshold_case_kwargs, mechanism_case_kwargs
export HydrodynamicConfig, biot_from_reynolds, detachment_coefficient, logspace
export run_nutrient_transfer_campaign, run_regime_map_campaign
export run_reynolds_threshold_campaign, run_reynolds_detachment_campaign

end
