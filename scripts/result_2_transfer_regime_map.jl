# Regime map sweep for the paper growth-extension model:
# computes R = Phi(T) / Phi(t_start) over (Bi, A0)

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "reduced_1d", "main_1d_paper.jl"))

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

function regime_label(R_post::Real; tol::Real=0.02)
    if R_post > 1 + tol
        return "growth"
    elseif R_post < 1 - tol
        return "death"
    else
        return "balanced"
    end
end

const CFG = (
    # geometry / numerics
    L0tilde = 0.8,
    N = 320,
    Tfinal = 30.0,
    dt = 2.0e-4,
    save_dt = 0.05,
    t_start = 6.0,
    live_frac0 = 0.6,
    k_d = 0.0,

    # reduced-model groups used in the paper
    Da = 180.0,
    Da_a = 0.0,
    Lambda = 0.2,
    Ktilde = 2.857142857142857,
    beta_Ab = 0.2222222222222222,
    PiE0 = 21.0,
    kill_power = 1.0,
    gamma_kill = 0.8,

    # sweep values
    Bi_values = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0],
    A_values  = [0.0, 0.005, 0.01, 0.02, 0.041, 0.08, 0.16],
)

const SAVE_EVERY_LOCAL = max(1, round(Int, CFG.save_dt / CFG.dt))
const OUTDIR = joinpath(REPO_ROOT, "results", "raw", "result_2_regime_map")
mkpath(OUTDIR)

summary_file = joinpath(OUTDIR, "test2_summary.csv")
open(summary_file, "w") do io
    println(io,
        "tag,Bi,A0,Da,L0tilde,live_frac0,Tfinal,t_start," *
        "Phi0,Phi_start,Phi_final,Phi_ratio_total,Phi_ratio_post,regime,L_final"
    )
end

for A0 in CFG.A_values
    for Bi in CFG.Bi_values
        tag = @sprintf("Bi_%0.3g_A_%0.3g", Bi, A0)
        @printf("\nRunning %s ...\n", tag)

        _, _, summary = run_case(;
            outdir = OUTDIR,
            tag = tag,

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
            C_top_robin_Bi = Bi,
            C_top_robin_ref = 1.0,
            A_top_bc = :robin,
            A_top_robin_Bi = Bi,
            A_top = 0.0,
            A_step_time = CFG.t_start,
            A_step_value = A0,

            save_final_profiles = false,
        )

        t_hist = summary.t_history
        Phi_hist = summary.Phi_history
        @assert !isempty(t_hist) "Empty timeseries history. Check save_every."

        Phi0 = CFG.live_frac0 * CFG.L0tilde
        Phi_start = linear_interp_at(CFG.t_start, t_hist, Phi_hist)
        Phi_final = summary.Phi_final
        Phi_ratio_total = Phi_final / max(Phi0, 1e-12)
        Phi_ratio_post = Phi_final / max(Phi_start, 1e-12)
        regime = regime_label(Phi_ratio_post)

        open(summary_file, "a") do io
            println(io,
                string(
                    tag, ",",
                    Bi, ",",
                    A0, ",",
                    CFG.Da, ",",
                    CFG.L0tilde, ",",
                    CFG.live_frac0, ",",
                    CFG.Tfinal, ",",
                    CFG.t_start, ",",
                    Phi0, ",",
                    Phi_start, ",",
                    Phi_final, ",",
                    Phi_ratio_total, ",",
                    Phi_ratio_post, ",",
                    regime, ",",
                    summary.L_final
                )
            )
        end

        @printf("  R = Phi(T)/Phi(t_start) = %.4f | regime = %s | L(T)=%.4f\n",
            Phi_ratio_post, regime, summary.L_final)
    end
end

println("\nDone.")
println("Summary saved in:")
println(summary_file)
