const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "reduced_1d", "main_1d_paper.jl"))

using Dates
using Printf

const CFG = (
    run_name = "re_da_threshold_parameter_sets_fixed_thickness",
    keep_case_files = false,

    # Fixed-thickness setting used for Result 3.
    fixed_thickness = true,
    L0tilde = 0.8,
    Lmin = 1.0e-8,
    N = 320,
    Tfinal = 30.0,
    dt = 2.0e-4,
    save_every = 500,
    tstart = 6.0,
    live_frac0 = 0.6,
    k_d = 0.0,

    # Parameters kept common across the two effective configurations.
    Lambda = 0.2,
    Da_a = 0.0,
    PiE0 = 21.0,
    beta_Ab = 0.2222222222222222,
    kill_power = 1.0,
    gamma_kill = 0.8,
    C_init = 1.0,
    A_init = 0.0,

    # Two physical/effective parameter configurations.
    # open_flow keeps the original high-Re parameterisation.
    # microfluidic_effective uses a low-Re laminar window and microchannel-like
    # transport scales. It is intended as a regime-test configuration, not as a
    # calibration to a specific microfluidic experiment.
    parameter_sets = (
        (
            name = :open_flow,
            label = "open-flow/high-Re baseline",
            Ktilde = 2.857142857142857,
            alpha_o = 1.1111111111111112,
            Sc_C = 1000.0,
            hb_over_lc = 1.25e-3,
            Da_values = [60.0, 180.0, 600.0],
            Re_values = [1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6],
            Re_transition_low = 3.0e4,
            Re_transition_high = 3.0e5,
        ),
        (
            name = :microfluidic_effective,
            label = "microfluidic-effective low-Re configuration",
            Ktilde = 2.0,
            alpha_o = 3.3333333333333335,
            Sc_C = 1000.0,
            hb_over_lc = 0.10,
            Da_values = [6.0, 18.0, 60.0],
            Re_values = [1.0e-1, 3.0e-1, 1.0e0, 3.0e0, 1.0e1, 3.0e1, 1.0e2, 3.0e2, 1.0e3],
            # Set the transition far above the explored range so that the
            # microfluidic window uses the laminar branch only.
            Re_transition_low = 1.0e99,
            Re_transition_high = 1.0e100,
        ),
    ),

    # Threshold search.
    bracket_low_factor = 0.08,
    bracket_high_factor = 4.0,
    bracket_points = 10,
    max_expand_steps = 4,
    expand_factor = 1.8,
    bisection_steps = 10,

    modes = (:nutrient_only, :fully_coupled),
)

gompertz_smoothstep(x::Float64) = x * x * (3.0 - 2.0 * x)

function smooth_transition_weight(Re::Float64, Re_lo::Float64, Re_hi::Float64)
    @assert Re_lo > 0.0 && Re_hi > Re_lo

    if Re <= Re_lo
        return 0.0
    elseif Re >= Re_hi
        return 1.0
    else
        ξ = (log10(Re) - log10(Re_lo)) / (log10(Re_hi) - log10(Re_lo))
        return gompertz_smoothstep(ξ)
    end
end

sherwood_laminar(Re::Float64, Sc::Float64) = 0.664 * Re^(1 / 2) * Sc^(1 / 3)
sherwood_turbulent(Re::Float64, Sc::Float64) = 0.037 * Re^(4 / 5) * Sc^(1 / 3)

function sherwood_effective(Re::Float64, Sc::Float64; Re_lo::Float64, Re_hi::Float64)
    Sh_lam = sherwood_laminar(Re, Sc)
    Sh_turb = sherwood_turbulent(Re, Sc)
    χ = smooth_transition_weight(Re, Re_lo, Re_hi)
    Sh_eff = (1.0 - χ) * Sh_lam + χ * Sh_turb

    return (
        Sh_eff = Sh_eff,
        Sh_lam = Sh_lam,
        Sh_turb = Sh_turb,
        chi = χ,
    )
end

function bi_from_re(
    Re::Float64;
    Sc::Float64,
    alpha::Float64,
    hb_over_lc::Float64,
    Re_lo::Float64,
    Re_hi::Float64,
)
    sh = sherwood_effective(Re, Sc; Re_lo = Re_lo, Re_hi = Re_hi)
    Bi = alpha * hb_over_lc * sh.Sh_eff
    return merge((Bi = Bi,), sh)
end

function bi_from_re_pset(Re::Float64, pset)
    return bi_from_re(
        Re;
        Sc = pset.Sc_C,
        alpha = pset.alpha_o,
        hb_over_lc = pset.hb_over_lc,
        Re_lo = pset.Re_transition_low,
        Re_hi = pset.Re_transition_high,
    )
end

function geom_range(lo::Float64, hi::Float64, n::Int)
    @assert lo > 0.0 && hi > lo && n >= 2
    return 10.0 .^ collect(range(log10(lo), log10(hi), length = n))
end

function linear_interp_at(x::Float64, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    n = length(xs)
    @assert n == length(ys) && n >= 1

    if x <= xs[1]
        return float(ys[1])
    elseif x >= xs[end]
        return float(ys[end])
    end

    idx = searchsortedfirst(xs, x)
    x2, y2 = float(xs[idx]), float(ys[idx])
    x1, y1 = float(xs[idx - 1]), float(ys[idx - 1])

    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end

function compute_R_post(summary, tstart::Float64)
    phi_start = if hasproperty(summary, :Phi_tstart) && isfinite(summary.Phi_tstart)
        summary.Phi_tstart
    else
        t_hist = collect(summary.t_history)
        phi_hist = collect(summary.Phi_history)
        @assert !isempty(t_hist) "Empty time history; increase save_every or inspect solver output."
        linear_interp_at(tstart, t_hist, phi_hist)
    end

    phi_final = summary.Phi_final
    R_post = phi_final / max(phi_start, 1e-14)

    return (
        phi_start = phi_start,
        phi_final = phi_final,
        R_post = R_post,
    )
end

function cleanup_case_files(outdir::String, tag::String)
    for fname in readdir(outdir)
        if startswith(fname, tag * "_")
            rm(joinpath(outdir, fname); force = true)
        end
    end
    return nothing
end

function evaluate_case(
    mode::Symbol,
    pset,
    Re::Float64,
    Da::Float64,
    A0::Float64;
    outdir::String,
    tag::String,
)
    map = bi_from_re_pset(Re, pset)
    Bi_C = map.Bi
    Bi_A = (mode === :fully_coupled) ? map.Bi : 0.0

    A_TOP_BC = if mode === :nutrient_only
        :dirichlet
    elseif mode === :fully_coupled
        :robin
    else
        error("Unknown mode: $mode")
    end

    _, _, summary = run_case(;
        outdir = outdir,
        tag = tag,
        save_final_profiles = false,
        fixed_thickness = CFG.fixed_thickness,

        L0tilde = CFG.L0tilde,
        Lmin = CFG.Lmin,
        N = CFG.N,
        Tfinal = CFG.Tfinal,
        dt = CFG.dt,
        save_every = CFG.save_every,

        Lambda = CFG.Lambda,
        Da = Da,
        Da_a = CFG.Da_a,
        beta_Ab = CFG.beta_Ab,
        Ktilde = pset.Ktilde,
        PiE0 = CFG.PiE0,
        kill_power = CFG.kill_power,
        gamma_kill = CFG.gamma_kill,
        live_frac0 = CFG.live_frac0,
        k_d = CFG.k_d,

        C_init = CFG.C_init,
        A_init = CFG.A_init,

        C_top_bc = :robin,
        C_top = 1.0,
        C_top_robin_Bi = Bi_C,
        C_top_robin_ref = 1.0,

        A_top_bc = A_TOP_BC,
        A_top = 0.0,
        A_top_robin_Bi = Bi_A,
        A_step_time = CFG.tstart,
        A_step_value = A0,
    )

    if CFG.fixed_thickness
        if abs(summary.L_final - CFG.L0tilde) > 1.0e-8
            error(
                "Fixed-thickness check failed: L_final=$(summary.L_final), " *
                "expected L0tilde=$(CFG.L0tilde). " *
                "Check main_1d_paper.jl fixed_thickness implementation."
            )
        end
    end

    resp = compute_R_post(summary, CFG.tstart)

    if !CFG.keep_case_files
        cleanup_case_files(outdir, tag)
    end

    return merge(map, resp, (
        scenario = String(pset.name),
        scenario_label = pset.label,
        mode = mode,
        Re = Re,
        Da = Da,
        A0 = A0,
        Ktilde = pset.Ktilde,
        alpha_o = pset.alpha_o,
        Sc_C = pset.Sc_C,
        hb_over_lc = pset.hb_over_lc,
        Bi_C = Bi_C,
        Bi_A = Bi_A,
        L_final = summary.L_final,
        Phi_final = summary.Phi_final,
    ))
end

function coarse_bracket(mode::Symbol, pset, Re::Float64, Da::Float64; outdir::String, case_id::Int)
    a0_ref = 1.0 / CFG.PiE0
    lo = CFG.bracket_low_factor * a0_ref
    hi = CFG.bracket_high_factor * a0_ref
    rows = NamedTuple[]

    for expand in 0:CFG.max_expand_steps
        Avals = vcat(0.0, geom_range(lo, hi, CFG.bracket_points))
        Rvals = Float64[]

        for (j, A0) in enumerate(Avals)
            tag = @sprintf("case_%03d_%s_coarse_e%d_j%d", case_id, String(pset.name), expand, j)
            ev = evaluate_case(mode, pset, Re, Da, A0; outdir = outdir, tag = tag)
            push!(Rvals, ev.R_post)

            push!(rows, (
                case_id = case_id,
                scenario = String(pset.name),
                scenario_label = pset.label,
                stage = "coarse",
                expand_step = expand,
                mode = String(mode),
                Re = Re,
                Da = Da,
                A0 = A0,
                R_post = ev.R_post,
                phi_start = ev.phi_start,
                phi_final = ev.phi_final,
                L_final = ev.L_final,
                Ktilde = ev.Ktilde,
                alpha_o = ev.alpha_o,
                Sc_C = ev.Sc_C,
                hb_over_lc = ev.hb_over_lc,
                Bi_C = ev.Bi_C,
                Bi_A = ev.Bi_A,
                Sh_eff = ev.Sh_eff,
                Sh_lam = ev.Sh_lam,
                Sh_turb = ev.Sh_turb,
                chi = ev.chi,
            ))
        end

        for j in eachindex(Avals)
            if abs(Rvals[j] - 1.0) < 1e-12
                return (status = "exact", Alo = Avals[j], Ahi = Avals[j], Rlo = Rvals[j], Rhi = Rvals[j], rows = rows)
            end
        end

        if Avals[1] == 0.0 && Rvals[1] < 1.0
            return (status = "already_below_zero", Alo = 0.0, Ahi = 0.0, Rlo = Rvals[1], Rhi = Rvals[1], rows = rows)
        end

        for j in 1:(length(Avals) - 1)
            R1, R2 = Rvals[j], Rvals[j + 1]
            if (R1 - 1.0) * (R2 - 1.0) < 0.0
                return (status = "bracketed", Alo = Avals[j], Ahi = Avals[j + 1], Rlo = R1, Rhi = R2, rows = rows)
            end
        end

        if all(r > 1.0 for r in Rvals)
            hi *= CFG.expand_factor
        elseif all(r < 1.0 for r in Rvals)
            lo /= CFG.expand_factor
        else
            return (status = "nonmonotone_or_unresolved", Alo = Avals[1], Ahi = Avals[end], Rlo = Rvals[1], Rhi = Rvals[end], rows = rows)
        end
    end

    return (status = "unresolved_after_expand", Alo = lo, Ahi = hi, Rlo = NaN, Rhi = NaN, rows = rows)
end

function bisection_threshold(
    mode::Symbol,
    pset,
    Re::Float64,
    Da::Float64,
    Alo::Float64,
    Ahi::Float64,
    Rlo::Float64,
    Rhi::Float64;
    outdir::String,
    case_id::Int,
)
    rows = NamedTuple[]
    lo, hi = Alo, Ahi
    rlo, rhi = Rlo, Rhi

    for k in 1:CFG.bisection_steps
        amid = 0.5 * (lo + hi)
        tag = @sprintf("case_%03d_%s_bisect_k%02d", case_id, String(pset.name), k)
        ev = evaluate_case(mode, pset, Re, Da, amid; outdir = outdir, tag = tag)
        rmid = ev.R_post

        push!(rows, (
            case_id = case_id,
            scenario = String(pset.name),
            scenario_label = pset.label,
            stage = "bisection",
            expand_step = k,
            mode = String(mode),
            Re = Re,
            Da = Da,
            A0 = amid,
            R_post = rmid,
            phi_start = ev.phi_start,
            phi_final = ev.phi_final,
            L_final = ev.L_final,
            Ktilde = ev.Ktilde,
            alpha_o = ev.alpha_o,
            Sc_C = ev.Sc_C,
            hb_over_lc = ev.hb_over_lc,
            Bi_C = ev.Bi_C,
            Bi_A = ev.Bi_A,
            Sh_eff = ev.Sh_eff,
            Sh_lam = ev.Sh_lam,
            Sh_turb = ev.Sh_turb,
            chi = ev.chi,
        ))

        if abs(rmid - 1.0) < 1e-12
            return (A0crit = amid, Rcrit = rmid, rows = rows)
        end

        if (rlo - 1.0) * (rmid - 1.0) <= 0.0
            hi = amid
            rhi = rmid
        else
            lo = amid
            rlo = rmid
        end
    end

    A0crit = 0.5 * (lo + hi)
    final_eval = evaluate_case(mode, pset, Re, Da, A0crit; outdir = outdir, tag = @sprintf("case_%03d_%s_final_recheck", case_id, String(pset.name)))

    return (A0crit = A0crit, Rcrit = final_eval.R_post, rows = rows)
end

function write_config(path::String)
    open(path, "w") do io
        println(io, "run_name = $(CFG.run_name)")
        println(io, "fixed_thickness = $(CFG.fixed_thickness)")
        println(io, "L0tilde = $(CFG.L0tilde)")
        println(io, "Lmin = $(CFG.Lmin)")
        println(io, "N = $(CFG.N)")
        println(io, "Tfinal = $(CFG.Tfinal)")
        println(io, "dt = $(CFG.dt)")
        println(io, "save_every = $(CFG.save_every)")
        println(io, "tstart = $(CFG.tstart)")
        println(io, "live_frac0 = $(CFG.live_frac0)")
        println(io, "k_d = $(CFG.k_d)")
        println(io, "Lambda = $(CFG.Lambda)")
        println(io, "Da_a = $(CFG.Da_a)")
        println(io, "PiE0 = $(CFG.PiE0)")
        println(io, "beta_Ab = $(CFG.beta_Ab)")
        println(io, "kill_power = $(CFG.kill_power)")
        println(io, "gamma_kill = $(CFG.gamma_kill)")
        println(io, "C_init = $(CFG.C_init)")
        println(io, "A_init = $(CFG.A_init)")
        println(io, "bracket_low_factor = $(CFG.bracket_low_factor)")
        println(io, "bracket_high_factor = $(CFG.bracket_high_factor)")
        println(io, "bracket_points = $(CFG.bracket_points)")
        println(io, "max_expand_steps = $(CFG.max_expand_steps)")
        println(io, "expand_factor = $(CFG.expand_factor)")
        println(io, "bisection_steps = $(CFG.bisection_steps)")
        println(io, "modes = $(collect(CFG.modes))")
        println(io, "")
        println(io, "parameter_sets:")
        for pset in CFG.parameter_sets
            println(io, "  name = $(pset.name)")
            println(io, "    label = $(pset.label)")
            println(io, "    Ktilde = $(pset.Ktilde)")
            println(io, "    alpha_o = $(pset.alpha_o)")
            println(io, "    Sc_C = $(pset.Sc_C)")
            println(io, "    hb_over_lc = $(pset.hb_over_lc)")
            println(io, "    Da_values = $(collect(pset.Da_values))")
            println(io, "    Re_values = $(collect(pset.Re_values))")
            println(io, "    Re_transition_low = $(pset.Re_transition_low)")
            println(io, "    Re_transition_high = $(pset.Re_transition_high)")
        end
    end
    return nothing
end

function append_eval_rows!(evals_csv::String, rows)
    open(evals_csv, "a") do io
        for row in rows
            println(
                io,
                "$(row.case_id),$(row.scenario),$(row.scenario_label),$(row.stage),$(row.expand_step),$(row.mode)," *
                "$(row.Re),$(row.Da),$(row.A0),$(row.R_post)," *
                "$(row.phi_start),$(row.phi_final),$(row.L_final)," *
                "$(row.Ktilde),$(row.alpha_o),$(row.Sc_C),$(row.hb_over_lc)," *
                "$(row.Bi_C),$(row.Bi_A),$(row.Sh_eff),$(row.Sh_lam),$(row.Sh_turb),$(row.chi)"
            )
        end
    end
    return nothing
end

function main()
    outdir = joinpath(REPO_ROOT, "results", "raw", "result_3_reynolds_threshold")
    if isdir(outdir)
        rm(outdir; recursive=true, force=true)
    end
    mkpath(outdir)

    surface_csv = joinpath(outdir, "re_da_threshold_surface.csv")
    evals_csv = joinpath(outdir, "re_da_threshold_evals.csv")
    mapping_csv = joinpath(outdir, "re_da_mapping.csv")
    config_txt = joinpath(outdir, "re_da_config.txt")

    cases_dir = joinpath(outdir, "cases")
    mkpath(cases_dir)
    write_config(config_txt)

    open(surface_csv, "w") do io
        println(
            io,
            "case_id,scenario,scenario_label,mode,fixed_thickness,Re,Da,Lambda,PiE0,Ktilde," *
            "alpha_o,Sc_C,hb_over_lc,Sh_eff,Sh_lam,Sh_turb,chi,Bi_C,Bi_A,status," *
            "A0crit,Alo,Ahi,Rlo,Rhi,phi_start_crit,phi_final_crit," *
            "Rcrit,L_final_crit,n_bisect_steps"
        )
    end

    open(evals_csv, "w") do io
        println(
            io,
            "case_id,scenario,scenario_label,stage,expand_step,mode,Re,Da,A0,R_post,phi_start," *
            "phi_final,L_final,Ktilde,alpha_o,Sc_C,hb_over_lc,Bi_C,Bi_A,Sh_eff,Sh_lam,Sh_turb,chi"
        )
    end

    open(mapping_csv, "w") do io
        println(io, "scenario,scenario_label,Re,Ktilde,alpha_o,Sc_C,hb_over_lc,Sh_lam,Sh_turb,chi,Sh_eff,Bi")
        for pset in CFG.parameter_sets, Re in pset.Re_values
            map = bi_from_re_pset(Re, pset)
            println(
                io,
                "$(String(pset.name)),$(pset.label),$(Re),$(pset.Ktilde),$(pset.alpha_o),$(pset.Sc_C)," *
                "$(pset.hb_over_lc),$(map.Sh_lam),$(map.Sh_turb),$(map.chi),$(map.Sh_eff),$(map.Bi)"
            )
        end
    end

    case_id = 0
    total_cases = sum(length(CFG.modes) * length(pset.Da_values) * length(pset.Re_values) for pset in CFG.parameter_sets)

    for pset in CFG.parameter_sets, mode in CFG.modes, Da in pset.Da_values, Re in pset.Re_values
        case_id += 1

        @printf(
            "\n[%d/%d] scenario=%s | mode=%s | Re=%.3g | Da=%.3g | Ktilde=%.3g\n",
            case_id,
            total_cases,
            String(pset.name),
            String(mode),
            Re,
            Da,
            pset.Ktilde,
        )

        coarse = coarse_bracket(mode, pset, Re, Da; outdir = cases_dir, case_id = case_id)
        append_eval_rows!(evals_csv, coarse.rows)

        status = coarse.status
        A0crit = NaN
        Rcrit = NaN
        phi_start_crit = NaN
        phi_final_crit = NaN
        L_final_crit = NaN
        n_bisect_steps = 0

        if status == "exact" || status == "already_below_zero"
            A0crit = coarse.Alo
            final_eval = evaluate_case(mode, pset, Re, Da, A0crit; outdir = cases_dir, tag = @sprintf("case_%03d_%s_finalcrit", case_id, String(pset.name)))
            phi_start_crit = final_eval.phi_start
            phi_final_crit = final_eval.phi_final
            Rcrit = final_eval.R_post
            L_final_crit = final_eval.L_final

        elseif status == "bracketed"
            bis = bisection_threshold(mode, pset, Re, Da, coarse.Alo, coarse.Ahi, coarse.Rlo, coarse.Rhi; outdir = cases_dir, case_id = case_id)
            A0crit = bis.A0crit
            Rcrit = bis.Rcrit
            n_bisect_steps = length(bis.rows)
            append_eval_rows!(evals_csv, bis.rows)

            final_eval = evaluate_case(mode, pset, Re, Da, A0crit; outdir = cases_dir, tag = @sprintf("case_%03d_%s_finalcrit", case_id, String(pset.name)))
            phi_start_crit = final_eval.phi_start
            phi_final_crit = final_eval.phi_final
            Rcrit = final_eval.R_post
            L_final_crit = final_eval.L_final
        end

        map = bi_from_re_pset(Re, pset)
        Bi_C = map.Bi
        Bi_A = (mode === :fully_coupled) ? map.Bi : 0.0

        open(surface_csv, "a") do io
            println(
                io,
                "$(case_id),$(String(pset.name)),$(pset.label),$(String(mode)),$(CFG.fixed_thickness),$(Re),$(Da),$(CFG.Lambda)," *
                "$(CFG.PiE0),$(pset.Ktilde),$(pset.alpha_o),$(pset.Sc_C),$(pset.hb_over_lc)," *
                "$(map.Sh_eff),$(map.Sh_lam),$(map.Sh_turb),$(map.chi)," *
                "$(Bi_C),$(Bi_A),$(status),$(A0crit),$(coarse.Alo)," *
                "$(coarse.Ahi),$(coarse.Rlo),$(coarse.Rhi),$(phi_start_crit)," *
                "$(phi_final_crit),$(Rcrit),$(L_final_crit),$(n_bisect_steps)"
            )
        end
    end

    @printf("\nSaved parameter-set Reynolds-Damköhler threshold sweep to: %s\n", outdir)
    @printf("  surface csv: %s\n", surface_csv)
    @printf("  evals csv  : %s\n", evals_csv)
    @printf("  mapping csv: %s\n", mapping_csv)
    @printf("  config txt : %s\n", config_txt)
    return nothing
end

main()
