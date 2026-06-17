#!/usr/bin/env julia

# Run the Julia scripts used to regenerate the computational results.
# Each result is executed in a separate Julia process so that scripts remain
# independent and reproducible from the repository root.

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function run_script(script_name::String)
    script_path = joinpath(@__DIR__, script_name)
    @info "Running $(script_name)"
    cmd = `$(Base.julia_cmd()) --project=$(REPO_ROOT) $(script_path)`
    run(cmd)
end

mkpath(joinpath(REPO_ROOT, "results"))

run_script("result_1_nutrient_transfer.jl")
run_script("result_1_asymptotic_profiles.jl")
run_script("result_2_transfer_regime_map.jl")
run_script("result_3_reynolds_threshold.jl")
run_script("result_4_stress_detachment.jl")

println("\nDone.")
println("Simulation outputs: ", joinpath(REPO_ROOT, "results", "raw"))
