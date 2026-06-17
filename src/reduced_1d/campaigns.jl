using DelimitedFiles
using Printf

function write_named_rows(path::String, header::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, join(row, ","))
        end
    end
    return path
end

function run_nutrient_transfer_campaign(;
    outdir::String="results/paper/nutrient_transfer",
    Bi_values::Vector{Float64}=[0.3, 1.0, 3.0, 10.0, 30.0, 100.0],
    kwargs...,
)
    rows = Any[]
    for BiC in Bi_values
        tag = @sprintf("nutrient_transfer_BiC_%0.4g", BiC)
        _, _, summary = run_case(;
            mechanism_case_kwargs(BiC=BiC)...,
            fixed_thickness=true,
            outdir=outdir,
            tag=replace(tag, "." => "p"),
            kwargs...,
        )
        push!(rows, (BiC, summary.R_final, summary.Phi_final, summary.Phi_tstart, summary.Cavg_history[end], summary.Aavg_history[end]))
    end
    write_named_rows(joinpath(outdir, "nutrient_transfer_summary.csv"),
                     "BiC,R_final,Phi_final,Phi_tstart,Cavg_final,Aavg_final", rows)
    return rows
end

function run_regime_map_campaign(;
    outdir::String="results/paper/regime_map",
    Bi_values::Vector{Float64}=logspace(0.1, 100.0, 31),
    A0_values::Vector{Float64}=collect(range(0.0, 0.12, length=31)),
    kwargs...,
)
    rows = Any[]
    for Bi in Bi_values, A0 in A0_values
        tag = @sprintf("regime_Bi_%0.4g_A_%0.4g", Bi, A0)
        _, _, summary = run_case(;
            threshold_case_kwargs(Bi=Bi, A0=A0)...,
            fixed_thickness=true,
            write_output=false,
            outdir=nothing,
            tag=replace(tag, "." => "p"),
            kwargs...,
        )
        push!(rows, (Bi, A0, summary.R_final))
    end
    mkpath(outdir)
    write_named_rows(joinpath(outdir, "regime_map.csv"), "Bi,A0,R_final", rows)
    return rows
end

function critical_A0_for_case(; target::Float64=1.0, A_low::Float64=0.0, A_high::Float64=0.20, tolerance::Float64=1.0e-3, maxiter::Int=30, kwargs...)
    base_kwargs = (; kwargs...)

    function f(A0)
        local_kwargs = merge(base_kwargs, (
            A_step_value=A0,
            write_output=false,
            outdir=nothing,
            fixed_thickness=true,
        ))
        _, _, summary = run_case(; local_kwargs...)
        return summary.R_final - target
    end

    lo = A_low
    hi = A_high
    flo = f(lo)
    fhi = f(hi)
    while sign(flo) == sign(fhi) && hi < 10.0
        hi *= 2.0
        fhi = f(hi)
    end
    if sign(flo) == sign(fhi)
        return NaN
    end

    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tolerance
            return mid
        elseif sign(fm) == sign(flo)
            lo = mid
            flo = fm
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

function run_reynolds_threshold_campaign(;
    outdir::String="results/paper/reynolds_threshold",
    configurations::Vector{Symbol}=[:open, :micro],
    nRe::Int=31,
    critical_tolerance::Float64=1.0e-3,
    critical_maxiter::Int=30,
    kwargs...,
)
    rows = Any[]
    for name in configurations
        config = HydrodynamicConfig(name)
        for Re in logspace(config.reynolds_min, config.reynolds_max, nRe)
            Bi = biot_from_reynolds(config, Re)

            Acrit_nutrient_only = critical_A0_for_case(;
                mechanism_case_kwargs(BiC=Bi)...,
                Ktilde=config.Ktilde,
                Da=config.Da_values[2],
                tolerance=critical_tolerance,
                maxiter=critical_maxiter,
                kwargs...,
            )

            Acrit_coupled = critical_A0_for_case(;
                threshold_case_kwargs(Bi=Bi, A0=0.04)...,
                Ktilde=config.Ktilde,
                Da=config.Da_values[2],
                tolerance=critical_tolerance,
                maxiter=critical_maxiter,
                kwargs...,
            )

            push!(rows, (String(name), Re, Bi, Acrit_nutrient_only, Acrit_coupled))
        end
    end
    mkpath(outdir)
    write_named_rows(joinpath(outdir, "reynolds_threshold.csv"),
                     "configuration,Re,Bi,Acrit_nutrient_only,Acrit_fully_coupled", rows)
    return rows
end

function run_reynolds_detachment_campaign(;
    outdir::String="results/paper/reynolds_detachment",
    configurations::Vector{Symbol}=[:open, :micro],
    nRe::Int=31,
    A0::Float64=0.041,
    kwargs...,
)
    rows = Any[]
    for name in configurations
        config = HydrodynamicConfig(name)
        σref = stress_reference(config)
        for Da_value in config.Da_values, Re in logspace(config.reynolds_min, config.reynolds_max, nRe)
            Bi = biot_from_reynolds(config, Re)
            kd_law = (t, L) -> detachment_coefficient(config, Re, L; σ_ref=σref)
            tag = replace(@sprintf("%s_Da_%0.4g_Re_%0.4g", String(name), Da_value, Re), "." => "p")

            _, _, summary = run_case(;
                C_top_bc=:robin,
                C_top_robin_Bi=Bi,
                C_top_robin_ref=1.0,
                A_top_bc=:dirichlet,
                A_step_time=6.0,
                A_step_value=A0,
                Ktilde=config.Ktilde,
                Da=Da_value,
                k_d_law=kd_law,
                fixed_thickness=false,
                save_final_profiles=false,
                outdir=joinpath(outdir, "timeseries"),
                tag=tag,
                kwargs...,
            )
            push!(rows, (String(name), Da_value, Re, Bi, summary.G_final, summary.R_final, summary.Phi_final, summary.L_final, summary.kd_history[end]))
        end
    end
    mkpath(outdir)
    write_named_rows(joinpath(outdir, "reynolds_detachment_summary.csv"),
                     "configuration,Da,Re,Bi,G_final,R_final,Phi_final,L_final,kd_final", rows)
    return rows
end
