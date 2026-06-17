const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "reduced_1d", "main_1d_paper.jl"))

using DelimitedFiles
using Printf

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

function first_crossing_time(t::AbstractVector{<:Real}, y::AbstractVector{<:Real}, level::Real)
    @assert length(t) == length(y)
    for i in eachindex(y)
        if y[i] <= level
            return float(t[i])
        end
    end
    return NaN
end

function write_timeseries_csv(filename::String, t_hist, t_start, Phi_hist, Phi_start, Cavg_hist, Aavg_hist, L_hist,
                              case_label::String, BiC::Real, A0::Real)
    post_window = max(last(t_hist) - t_start, 1e-12)
    open(filename, "w") do io
        println(io, "case,label,BiC,A0,t,t_post,tau_post,L,Phi,Phi_ratio_post,Cavg_bio,Aavg_bio")
        for i in eachindex(t_hist)
            t_now = float(t_hist[i])
            t_post = t_now - t_start
            tau_post = t_post / post_window
            phi_now = float(Phi_hist[i])
            phi_ratio_post = phi_now / max(Phi_start, 1e-12)
            cavg = float(Cavg_hist[i])
            aavg = float(Aavg_hist[i])
            Lnow = float(L_hist[i])
            println(io,
                string(case_label, ",", case_label, ",", BiC, ",", A0, ",", t_now, ",", t_post, ",",
                       tau_post, ",", Lnow, ",", phi_now, ",", phi_ratio_post, ",", cavg, ",", aavg)
            )
        end
    end
end

function write_profile_csv(filename::String, z, C, A, e, phi, v)
    open(filename, "w") do io
        println(io, "z,C,A,e,phi,v")
        data = hcat(z, C, A, e, phi, v)
        writedlm(io, data, ',')
    end
end

const CFG = (
    outdir = joinpath(REPO_ROOT, "results", "raw", "result_1_nutrient_transfer"),

    L0tilde = 0.8,
    N = 320,
    Tfinal = 30.0,
    dt = 2.0e-4,
    save_dt = 0.05,
    t_start = 6.0,
    snapshot_post_frac = 0.50,
    live_frac0 = 0.6,
    k_d = 0.0,

    Da = 180.0,
    Da_a = 0.0,
    Lambda = 0.2,
    Ktilde = 2.857142857142857,
    beta_Ab = 0.2222222222222222,
    PiE0 = 54.0,
    kill_power = 1.0,
    gamma_kill = 0.8,

    A0_fixed = 0.11,
    antibiotic_bc_mode = :dirichlet,
    BiA_fixed = 30.0,
    C_bulk_ref = 1.0,

    cases = [
        (key = "low",  label = "Low nutrient transfer",          BiC = 0.1),
        (key = "mid",  label = "Intermediate nutrient transfer", BiC = 1.0),
        (key = "high", label = "High nutrient transfer",         BiC = 30.0),
    ],
)

const SAVE_EVERY_LOCAL = max(1, round(Int, CFG.save_dt / CFG.dt))
const T_SNAPSHOT = CFG.t_start + CFG.snapshot_post_frac * (CFG.Tfinal - CFG.t_start)
mkpath(CFG.outdir)

summary_file = joinpath(CFG.outdir, "mechanism_summary.csv")
open(summary_file, "w") do io
    println(io,
        "case,label,BiC,A0,antibiotic_bc,BiA,L0tilde,Da,live_frac0,Tfinal," *
        "t_start,t_snapshot,Phi0,Phi_start,Phi_snapshot,Phi_final," *
        "Phi_ratio_snapshot,Phi_ratio_post,t_half_post,L_snapshot,L_final"
    )
end

for case in CFG.cases
    tag_full = "mechanism_$(case.key)_full"
    @printf("\nRunning %s (Bi_C = %.3g) ...\n", tag_full, case.BiC)

    A_TOP_BC = if CFG.antibiotic_bc_mode === :dirichlet
        :dirichlet
    elseif CFG.antibiotic_bc_mode === :robin
        :robin
    else
        error("Unknown antibiotic_bc_mode")
    end
    A_TOP_ROBIN_BI = CFG.antibiotic_bc_mode === :robin ? CFG.BiA_fixed : 0.0

    # Full run to Tfinal for timeseries and final state
    _, _, summary = run_case(;
        outdir = CFG.outdir,
        tag = tag_full,
        L0tilde = CFG.L0tilde,
        N = CFG.N,
        Tfinal = CFG.Tfinal,
        dt = CFG.dt,
        save_every = SAVE_EVERY_LOCAL,
        Lambda = CFG.Lambda,
        Da = CFG.Da,
        Da_a = CFG.Da_a,
        beta_Ab = CFG.beta_Ab,
        Ktilde = CFG.Ktilde,
        PiE0 = CFG.PiE0,
        kill_power = CFG.kill_power,
        gamma_kill = CFG.gamma_kill,
        live_frac0 = CFG.live_frac0,
        k_d = CFG.k_d,
        fixed_thickness = true,
        C_top_bc = :robin,
        C_top_robin_Bi = case.BiC,
        C_top_robin_ref = CFG.C_bulk_ref,
        A_top_bc = A_TOP_BC,
        A_top = 0.0,
        A_top_robin_Bi = A_TOP_ROBIN_BI,
        A_step_time = CFG.t_start,
        A_step_value = CFG.A0_fixed,
    )

    t_hist = summary.t_history
    Phi_hist = summary.Phi_history
    Cavg_hist = summary.Cavg_history
    Aavg_hist = summary.Aavg_history
    L_hist = summary.L_history
    @assert !isempty(t_hist)

    Phi0 = CFG.live_frac0 * CFG.L0tilde
    Phi_start = linear_interp_at(CFG.t_start, t_hist, Phi_hist)
    Phi_final = summary.Phi_final
    Phi_ratio_post = Phi_final / max(Phi_start, 1e-12)

    # Snapshot run terminating at T_snapshot to recover profiles there
    tag_snap = "mechanism_$(case.key)_snapshot"
    _, _, summary_snap = run_case(;
        outdir = CFG.outdir,
        tag = tag_snap,
        L0tilde = CFG.L0tilde,
        N = CFG.N,
        Tfinal = T_SNAPSHOT,
        dt = CFG.dt,
        save_every = SAVE_EVERY_LOCAL,
        Lambda = CFG.Lambda,
        Da = CFG.Da,
        Da_a = CFG.Da_a,
        beta_Ab = CFG.beta_Ab,
        Ktilde = CFG.Ktilde,
        PiE0 = CFG.PiE0,
        kill_power = CFG.kill_power,
        gamma_kill = CFG.gamma_kill,
        live_frac0 = CFG.live_frac0,
        k_d = CFG.k_d,
        fixed_thickness = true,
        C_top_bc = :robin,
        C_top_robin_Bi = case.BiC,
        C_top_robin_ref = CFG.C_bulk_ref,
        A_top_bc = A_TOP_BC,
        A_top = 0.0,
        A_top_robin_Bi = A_TOP_ROBIN_BI,
        A_step_time = CFG.t_start,
        A_step_value = CFG.A0_fixed,
    )

    Phi_snapshot = summary_snap.Phi_final
    Phi_ratio_snapshot = Phi_snapshot / max(Phi_start, 1e-12)
    t_half_post = first_crossing_time(t_hist, Phi_hist ./ max(Phi_start, 1e-12), 0.5)

    tsfile = joinpath(CFG.outdir, "mechanism_$(case.key)_timeseries.csv")
    write_timeseries_csv(tsfile, t_hist, CFG.t_start, Phi_hist, Phi_start, Cavg_hist, Aavg_hist, L_hist,
                         case.label, case.BiC, CFG.A0_fixed)

    snapfile = joinpath(CFG.outdir, "mechanism_$(case.key)_snapshot_profiles.csv")
    write_profile_csv(snapfile, summary_snap.z, summary_snap.C_final, summary_snap.A_final,
                      summary_snap.e_final, summary_snap.phi_final, summary_snap.v_final)

    open(summary_file, "a") do io
        println(io,
            string(case.key, ",", case.label, ",", case.BiC, ",", CFG.A0_fixed, ",",
                   CFG.antibiotic_bc_mode, ",", (CFG.antibiotic_bc_mode === :robin ? CFG.BiA_fixed : 0.0), ",",
                   CFG.L0tilde, ",", CFG.Da, ",", CFG.live_frac0, ",", CFG.Tfinal, ",",
                   CFG.t_start, ",", T_SNAPSHOT, ",", Phi0, ",", Phi_start, ",", Phi_snapshot, ",", Phi_final, ",",
                   Phi_ratio_snapshot, ",", Phi_ratio_post, ",", t_half_post, ",", summary_snap.L_final, ",", summary.L_final)
        )
    end

    @printf("  R = Phi(T)/Phi(t_start) = %.4f | L(T)=%.4f | snapshot at t=%.2f\n",
            Phi_ratio_post, summary.L_final, T_SNAPSHOT)
end

println("\nDone.")
println("Summary:    " * summary_file)
println("Output dir: " * CFG.outdir)
