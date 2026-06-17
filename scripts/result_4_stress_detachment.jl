using Dates
using Printf
using DelimitedFiles

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "src", "reduced_1d", "main_1d_paper.jl"))

using .Grid1D
using .Transport1D
using .ConstraintNu1D

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
const CFG = (
    run_name = "re_stress_coupled_parameter_sets_v2",

    # Common numerical/protocol settings.
    L0tilde = 0.8,
    N = 320,
    Tfinal = 30.0,
    dt = 2.0e-4,
    save_every = 1000,
    tstart = 6.0,
    live_frac0 = 0.6,
    Lmin = 1.0e-4,

    # Common biological/treatment parameters.
    Lambda = 0.2,
    beta_Ab = 0.2222222222222222,
    Da_a = 0.0,
    PiE0 = 21.0,
    kill_power = 1.0,
    gamma_kill = 0.8,
    A0 = 0.041,

    # Stress-based detachment law.
    kdmax = 0.013,
    q_stress = 0.5,

    # Normalization strategy.
    # :scenario_reference_clogging means:
    #   sigma_eff_ref = max_Re sigma0(Re) * amplification(phi_c_ref_norm)
    stress_normalization = :scenario_reference_clogging,

    # Two physical/effective parameter configurations.
    parameter_sets = (
        (
            name = :open_flow,
            label = "open-flow/high-Re baseline",

            # Biofilm-domain/reaction-transfer parameters.
            Ktilde = 2.857142857142857,
            alpha_o = 1.1111111111111112,
            Sc_C = 1000.0,
            hb_over_lc = 1.25e-3,
            Da_values = [60.0, 180.0, 600.0],
            Re_values = [
                1.0e3,
                3.0e3,
                1.0e4,
                3.0e4,
                6.0e4,
                1.0e5,
                2.0e5,
                3.0e5,
                5.0e5,
                7.0e5,
                1.0e6,
            ],
            Re_transition_low = 3.0e4,
            Re_transition_high = 3.0e5,

            # Stress proxy.
            stress_model = :open_wall_shear,
            rho = 1000.0,
            nu = 1.0e-6,
            Lc_stress_m = 1.0e-2,
            H_channel_m = NaN,
            C_square = NaN,

            # No clogging amplification for open-flow/high-Re configuration.
            clogging_model = :none,
            chi_H = 0.0,
            phi_c_max = 0.0,
            phi_c_ref_norm = 0.0,
        ),
        (
            name = :microfluidic_effective,
            label = "microfluidic-effective low-Re configuration",

            # Biofilm-domain/reaction-transfer parameters.
            Ktilde = 2.0,
            alpha_o = 3.3333333333333335,
            Sc_C = 1000.0,
            hb_over_lc = 0.10,
            Da_values = [6.0, 18.0, 60.0],
            Re_values = [
                1.0e-1,
                3.0e-1,
                1.0e0,
                3.0e0,
                1.0e1,
                3.0e1,
                1.0e2,
                3.0e2,
                1.0e3,
            ],
            Re_transition_low = 1.0e99,
            Re_transition_high = 1.0e100,

            # Confined square-channel stress proxy:
            #   sigma0(Re) ~= C_square * rho * nu^2 / H^2 * Re
            # with C_square ~= 7.1 for a clean square channel.
            stress_model = :confined_square_channel,
            rho = 1000.0,
            nu = 1.0e-6,
            Lc_stress_m = NaN,
            H_channel_m = 100.0e-6,
            C_square = 7.1,

            # Effective clogging fraction from 1D thickness:
            #   phi_c(Ltilde) = min(chi_H * Ltilde, phi_c_max)
            #
            # chi_H controls how strongly the 1D thickness maps into channel
            # obstruction. For example, chi_H = 0.10 gives phi_c = 0.10 when
            # Ltilde = 1.
            #
            # phi_c_max is only a numerical cap.
            # phi_c_ref_norm is the reference clogging level used to normalize
            # the stress-based detachment law.
            clogging_model = :linear_from_thickness,
            chi_H = 0.10,
            phi_c_max = 0.95,
            phi_c_ref_norm = 0.20,
        ),
    ),

    # Final setting used for Result 4: nutrient transfer depends on Re, antibiotic
    # is imposed at the biofilm surface to isolate transport-detachment competition.
    mode = :nutrient_only,
)

# ------------------------------------------------------------------------------
# Transport mapping helpers
# ------------------------------------------------------------------------------
smoothstep(x::Float64) = x * x * (3.0 - 2.0 * x)

function transition_weight(Re::Float64, Re_lo::Float64, Re_hi::Float64)
    if Re <= Re_lo
        return 0.0
    elseif Re >= Re_hi
        return 1.0
    else
        xi = (log10(Re) - log10(Re_lo)) / (log10(Re_hi) - log10(Re_lo))
        return smoothstep(xi)
    end
end

sherwood_laminar(Re::Float64, Sc::Float64) = 0.664 * Re^(1/2) * Sc^(1/3)
sherwood_turbulent(Re::Float64, Sc::Float64) = 0.037 * Re^(4/5) * Sc^(1/3)

function sherwood_effective(Re::Float64, Sc::Float64; Re_lo::Float64, Re_hi::Float64)
    Sh_lam = sherwood_laminar(Re, Sc)
    Sh_turb = sherwood_turbulent(Re, Sc)
    chi = transition_weight(Re, Re_lo, Re_hi)
    Sh_eff = (1.0 - chi) * Sh_lam + chi * Sh_turb
    return (Sh_eff = Sh_eff, Sh_lam = Sh_lam, Sh_turb = Sh_turb, chi = chi)
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

# ------------------------------------------------------------------------------
# Stress-based detachment helpers
# ------------------------------------------------------------------------------
function skin_friction_flat_plate(Re::Float64; Re_lo::Float64, Re_hi::Float64)
    Re_safe = max(Re, eps(Float64))
    Cf_lam = 1.328 / sqrt(Re_safe)
    Cf_turb = 0.074 / Re_safe^(1/5)
    chi = transition_weight(Re_safe, Re_lo, Re_hi)
    Cf = (1.0 - chi) * Cf_lam + chi * Cf_turb
    return (Cf = Cf, Cf_lam = Cf_lam, Cf_turb = Cf_turb, chi_stress = chi)
end

function clean_stress_pset(Re::Float64, pset)
    if pset.stress_model === :confined_square_channel
        # Clean-channel stress for a square microchannel.
        sigma0 = pset.C_square * pset.rho * pset.nu^2 / pset.H_channel_m^2 * Re

        return (
            sigma0 = sigma0,
            U_ref = Re * pset.nu / pset.H_channel_m,
            stress_aux_1 = pset.C_square,
            stress_aux_2 = NaN,
            stress_model = String(pset.stress_model),
        )

    elseif pset.stress_model === :open_wall_shear
        # Wall-shear proxy:
        #   tau_w = 0.5 rho U^2 Cf,
        #   U = Re nu / Lc.
        cf = skin_friction_flat_plate(
            Re;
            Re_lo = pset.Re_transition_low,
            Re_hi = pset.Re_transition_high,
        )

        U = Re * pset.nu / pset.Lc_stress_m
        tau = 0.5 * pset.rho * U^2 * cf.Cf

        return (
            sigma0 = tau,
            U_ref = U,
            stress_aux_1 = cf.Cf,
            stress_aux_2 = cf.chi_stress,
            stress_model = String(pset.stress_model),
        )

    else
        error("Unknown stress_model: $(pset.stress_model)")
    end
end

function effective_clogging_fraction(Ltilde::Float64, pset)
    if pset.clogging_model === :none
        return 0.0

    elseif pset.clogging_model === :linear_from_thickness
        return clamp(pset.chi_H * Ltilde, 0.0, pset.phi_c_max)

    elseif pset.clogging_model === :square_wall_layer
        # Optional alternative for a wall layer in a square channel:
        #   free area ratio = (1 - 2 b/H)^2
        #   phi_c = 1 - free area ratio
        b_over_H = clamp(pset.chi_H * Ltilde, 0.0, 0.499999)
        phi_c = 1.0 - (1.0 - 2.0 * b_over_H)^2
        return clamp(phi_c, 0.0, pset.phi_c_max)

    else
        error("Unknown clogging_model: $(pset.clogging_model)")
    end
end

function clogging_amplification_from_phi(phi_c::Float64)
    phi_safe = clamp(phi_c, 0.0, 0.999999)
    return (1.0 - phi_safe)^(-2)
end

function stress_effective_pset(Re::Float64, Ltilde::Float64, pset)
    base = clean_stress_pset(Re, pset)
    phi_c = effective_clogging_fraction(Ltilde, pset)

    amp = if pset.clogging_model === :none
        1.0
    else
        clogging_amplification_from_phi(phi_c)
    end

    sigma_eff = base.sigma0 * amp

    return merge(base, (
        phi_c = phi_c,
        clogging_amplification = amp,
        sigma_eff = sigma_eff,
    ))
end

function sigma_eff_norm_pset(pset)
    # Corrected normalization.
    #
    # Old version:
    #   sigma_norm = max_Re sigma0(Re) * (1 - phi_c_max)^(-2)
    #
    # This made detachment too weak when phi_c_max was only a numerical cap.
    #
    # New version:
    #   sigma_norm = max_Re sigma0(Re) * (1 - phi_c_ref_norm)^(-2)
    #
    # phi_c_ref_norm is a physically relevant reference clogging level.
    vals = Float64[]

    for Re in pset.Re_values
        base = clean_stress_pset(Re, pset)

        amp_ref = if pset.clogging_model === :none
            1.0
        else
            @assert 0.0 <= pset.phi_c_ref_norm < 1.0
            clogging_amplification_from_phi(pset.phi_c_ref_norm)
        end

        push!(vals, base.sigma0 * amp_ref)
    end

    return maximum(vals)
end

function sigma_eff_ref_possible_pset(Re::Float64, pset)
    # Diagnostic only: effective stress at the reference normalization clogging.
    base = clean_stress_pset(Re, pset)

    amp_ref = if pset.clogging_model === :none
        1.0
    else
        clogging_amplification_from_phi(pset.phi_c_ref_norm)
    end

    return base.sigma0 * amp_ref
end

function sigma_eff_cap_possible_pset(Re::Float64, pset)
    # Diagnostic only: effective stress at the numerical cap.
    base = clean_stress_pset(Re, pset)

    amp_cap = if pset.clogging_model === :none
        1.0
    else
        clogging_amplification_from_phi(pset.phi_c_max)
    end

    return base.sigma0 * amp_cap
end

function kd_from_effective_stress(
    sigma_eff::Float64,
    sigma_norm::Float64,
    kdmax::Float64,
    q::Float64,
)
    @assert sigma_norm > 0.0

    # The clamp prevents k_d from exceeding kdmax. This keeps kdmax interpretable
    # as the maximum reduced detachment strength in the explored scenario.
    ratio = clamp(sigma_eff / sigma_norm, 0.0, 1.0)
    return kdmax * ratio^q
end

function make_kd_callback(Re::Float64, pset, sigma_norm::Float64)
    return function (Ltilde::Float64, t::Float64)
        sd = stress_effective_pset(Re, Ltilde, pset)
        kd = kd_from_effective_stress(sd.sigma_eff, sigma_norm, CFG.kdmax, CFG.q_stress)

        return merge(sd, (
            k_d = kd,
            sigma_eff_norm = sigma_norm,
            q_stress = CFG.q_stress,
            kdmax = CFG.kdmax,
        ))
    end
end

# ------------------------------------------------------------------------------
# Local dynamic-kd solver
# ------------------------------------------------------------------------------
function run_case_stress_detachment(;
    L0tilde::Float64 = CFG.L0tilde,
    N::Int = CFG.N,
    Tfinal::Float64 = CFG.Tfinal,
    dt::Float64 = CFG.dt,
    save_every::Int = CFG.save_every,
    Lambda::Float64 = CFG.Lambda,
    Da::Float64,
    Da_a::Float64 = CFG.Da_a,
    beta_Ab::Float64 = CFG.beta_Ab,
    Ktilde::Float64,
    PiE0::Float64 = CFG.PiE0,
    kill_power::Float64 = CFG.kill_power,
    gamma_kill::Float64 = CFG.gamma_kill,
    live_frac0::Float64 = CFG.live_frac0,
    Lmin::Float64 = CFG.Lmin,
    C_init::Float64 = 0.0,
    A_init::Float64 = 0.0,
    C_top_bc::Symbol = :robin,
    C_top::Float64 = 1.0,
    C_top_robin_Bi::Float64,
    C_top_robin_ref::Float64 = 1.0,
    A_top_bc::Symbol = :dirichlet,
    A_top::Float64 = 0.0,
    A_top_robin_Bi::Float64 = C_top_robin_Bi,
    A_step_time::Float64 = CFG.tstart,
    A_step_value::Float64 = CFG.A0,
    outdir::String,
    tag::String,
    save_final_profiles::Bool = false,
    fixed_thickness::Bool = false,
    k_d_callback,
)
    @assert C_top_bc in (:dirichlet, :robin)
    @assert A_top_bc in (:dirichlet, :robin)
    @assert 0.0 < live_frac0 <= 1.0

    L = L0tilde
    grid = Grid1D.make_grid(L, N)
    z = grid.z
    dz = grid.dz

    C = fill(C_init, N)
    A = fill(A_init, N)
    phi = fill(live_frac0, N)

    e = zeros(N)
    v = zeros(N)
    R_C = zeros(N)
    R_A = zeros(N)
    R_phi = zeros(N)
    D_C = fill(1.0, N)
    D_A = fill(beta_Ab, N)
    D_zero = fill(0.0, N)
    v_zero = fill(0.0, N)

    t_hist = Float64[]
    L_hist = Float64[]
    Phi_hist = Float64[]
    vtop_hist = Float64[]
    B_hist = Float64[]
    Cavg_hist = Float64[]
    Aavg_hist = Float64[]
    kd_hist = Float64[]
    sigma0_hist = Float64[]
    sigmaeff_hist = Float64[]
    phic_hist = Float64[]
    clogamp_hist = Float64[]

    Phi_tstart = NaN
    t = 0.0
    istep = 0

    mkpath(outdir)

    while t < Tfinal - 1e-15
        istep += 1
        t += dt

        Transport1D.monod!(e, C, Ktilde)

        if fixed_thickness
            fill!(v, 0.0)
        else
            ConstraintNu1D.compute_v!(v, e, phi, Lambda, dz)
        end

        @inbounds for i in 1:N
            R_C[i] = -Da * e[i] * phi[i]
            R_A[i] = -Da_a * A[i] * phi[i]
            R_phi[i] =
                Lambda * e[i] * phi[i] -
                Lambda * PiE0 * (e[i]^kill_power) * A[i] * (phi[i]^gamma_kill)
        end

        # Solutes on current physical grid.
        if C_top_bc === :robin
            Transport1D.step_scalar_imex!(C, v_zero, D_C, dt, dz;
                reaction_rhs = R_C,
                top_dirichlet = false,
                top_robin = true,
                top_robin_coeff = C_top_robin_Bi,
                top_robin_ref = C_top_robin_ref,
            )
        else
            Transport1D.step_scalar_imex!(C, v_zero, D_C, dt, dz;
                reaction_rhs = R_C,
                top_dirichlet = true,
                top_value = C_top,
            )
        end

        A_bc_now = (t >= A_step_time) ? A_step_value : A_top

        if A_top_bc === :robin
            Transport1D.step_scalar_imex!(A, v_zero, D_A, dt, dz;
                reaction_rhs = R_A,
                top_dirichlet = false,
                top_robin = true,
                top_robin_coeff = A_top_robin_Bi,
                top_robin_ref = A_bc_now,
            )
        else
            Transport1D.step_scalar_imex!(A, v_zero, D_A, dt, dz;
                reaction_rhs = R_A,
                top_dirichlet = true,
                top_value = A_bc_now,
            )
        end

        # Biomass with growth-induced advection.
        Transport1D.step_scalar_imex!(phi, v, D_zero, dt, dz;
            reaction_rhs = R_phi,
            top_dirichlet = false,
            advect_top_outflow = true,
        )

        @inbounds for i in 1:N
            C[i] = max(C[i], 0.0)
            A[i] = max(A[i], 0.0)
            phi[i] = max(phi[i], 0.0)
        end

        kd_diag = k_d_callback(L, t)
        kd_now = fixed_thickness ? 0.0 : kd_diag.k_d

        # Update moving boundary only when the moving-domain model is active.
        if !fixed_thickness
            Ldot = v[end] - kd_now * L^2
            Lnew = max(L + dt * Ldot, Lmin)

            if abs(Lnew - L) > 1e-14
                z_old = z
                C_old = copy(C)
                A_old = copy(A)
                phi_old = copy(phi)

                grid = Grid1D.make_grid(Lnew, N)
                z = grid.z
                dz = grid.dz

                C = linear_remap(z_old, C_old, z)
                A = linear_remap(z_old, A_old, z)
                phi = linear_remap(z_old, phi_old, z)

                e = similar(C)
                v = zeros(N)
                R_C = zeros(N)
                R_A = zeros(N)
                R_phi = zeros(N)
                D_C = fill(1.0, N)
                D_A = fill(beta_Ab, N)
                D_zero = fill(0.0, N)
                v_zero = fill(0.0, N)
                L = Lnew
            else
                L = Lnew
            end
        else
            L = L0tilde
        end

        if isnan(Phi_tstart) && t >= A_step_time
            Phi_tstart = Grid1D.trapz_uniform(phi, dz)
        end

        if istep % save_every == 0
            Phi = Grid1D.trapz_uniform(phi, dz)
            Cavg = Grid1D.trapz_uniform(C, dz) / L
            Aavg = Grid1D.trapz_uniform(A, dz) / L
            B = isnan(Phi_tstart) ? NaN : (Phi / max(Phi_tstart, 1e-12))
            kd_save = k_d_callback(L, t)

            push!(t_hist, t)
            push!(L_hist, L)
            push!(Phi_hist, Phi)
            push!(vtop_hist, v[end])
            push!(B_hist, B)
            push!(Cavg_hist, Cavg)
            push!(Aavg_hist, Aavg)
            push!(kd_hist, kd_save.k_d)
            push!(sigma0_hist, kd_save.sigma0)
            push!(sigmaeff_hist, kd_save.sigma_eff)
            push!(phic_hist, kd_save.phi_c)
            push!(clogamp_hist, kd_save.clogging_amplification)
        end
    end

    Phi_final = Grid1D.trapz_uniform(phi, dz)
    R_final = isnan(Phi_tstart) ? NaN : (Phi_final / max(Phi_tstart, 1e-12))
    final_diag = k_d_callback(L, t)

    summary = (
        L_final = L,
        Phi_final = Phi_final,
        Phi_tstart = Phi_tstart,
        R_final = R_final,
        k_d_final = final_diag.k_d,
        sigma0_final = final_diag.sigma0,
        sigma_eff_final = final_diag.sigma_eff,
        sigma_eff_norm = final_diag.sigma_eff_norm,
        phi_c_final = final_diag.phi_c,
        clogging_amplification_final = final_diag.clogging_amplification,
        t_history = t_hist,
        L_history = L_hist,
        Phi_history = Phi_hist,
        vtop_history = vtop_hist,
        B_history = B_hist,
        Cavg_history = Cavg_hist,
        Aavg_history = Aavg_hist,
        kd_history = kd_hist,
        sigma0_history = sigma0_hist,
        sigmaeff_history = sigmaeff_hist,
        phic_history = phic_hist,
        clogamp_history = clogamp_hist,
        z = z,
        C_final = copy(C),
        A_final = copy(A),
        phi_final = copy(phi),
        v_final = copy(v),
        e_final = copy(e),
    )

    tsfile = joinpath(outdir, "$(tag)_timeseries.csv")
    open(tsfile, "w") do io
        println(io, "t,L,Phi,v_top,B,Cavg,Aavg,k_d,sigma0,sigma_eff,phi_c,clogging_amplification")
        for i in eachindex(t_hist)
            println(
                io,
                "$(t_hist[i]),$(L_hist[i]),$(Phi_hist[i]),$(vtop_hist[i]),$(B_hist[i]),$(Cavg_hist[i]),$(Aavg_hist[i]),$(kd_hist[i]),$(sigma0_hist[i]),$(sigmaeff_hist[i]),$(phic_hist[i]),$(clogamp_hist[i])",
            )
        end
    end

    if save_final_profiles
        prof = hcat(z, C, A, phi, v, e)
        writedlm(joinpath(outdir, "$(tag)_final_profiles.csv"), prof, ',')
    end

    return z, L, summary
end

# ------------------------------------------------------------------------------
# Output helpers
# ------------------------------------------------------------------------------
append_csv_row(io, vals) = println(io, join(string.(vals), ","))

function classify_trend(vals::Vector{Float64}; atol::Float64 = 1e-8)
    if length(vals) <= 1
        return "degenerate"
    end

    nondecreasing = true
    nonincreasing = true

    for i in 2:length(vals)
        nondecreasing &= (vals[i] >= vals[i - 1] - atol)
        nonincreasing &= (vals[i] <= vals[i - 1] + atol)
    end

    if nondecreasing
        return "monotone_nondecreasing"
    elseif nonincreasing
        return "monotone_nonincreasing"
    end

    imax = argmax(vals)
    left_ok = true
    right_ok = true

    for i in 2:imax
        left_ok &= (vals[i] >= vals[i - 1] - atol)
    end

    for i in imax+1:length(vals)
        right_ok &= (vals[i] <= vals[i - 1] + atol)
    end

    if left_ok && right_ok && imax != 1 && imax != length(vals)
        return "single_peak_nonmonotone"
    end

    return "other_nonmonotone"
end

function peak_summary(Re_vals::Vector{Float64}, vals::Vector{Float64})
    i = argmax(vals)
    return (peak_idx = i, peak_Re = Re_vals[i], peak_val = vals[i])
end

safe_minimum(v::Vector{Float64}) = isempty(v) ? NaN : minimum(v)
safe_maximum(v::Vector{Float64}) = isempty(v) ? NaN : maximum(v)
safe_mean(v::Vector{Float64}) = isempty(v) ? NaN : sum(v) / length(v)

# ------------------------------------------------------------------------------
# Output folder and CSV headers
# ------------------------------------------------------------------------------
OUTDIR = joinpath(REPO_ROOT, "results", "raw", "result_4_stress_detachment")
if isdir(OUTDIR)
    rm(OUTDIR; recursive=true, force=true)
end
mkpath(OUTDIR)

summary_csv = joinpath(OUTDIR, "re_stress_coupled_parameter_sets_summary.csv")
figure_csv  = joinpath(OUTDIR, "re_stress_coupled_parameter_sets_figure_ready.csv")
mapping_csv = joinpath(OUTDIR, "re_stress_coupled_parameter_sets_mapping.csv")

open(summary_csv, "w") do io
    println(
        io,
        "scenario,scenario_label,mode,Da,Re,Ktilde,alpha_o,Sc_C,hb_over_lc,BiC,Sh_eff,Sh_lam,Sh_turb,chi,stress_model,clogging_model,chi_H,phi_c_max,phi_c_ref_norm,sigma0_clean,sigma_eff_norm,kdmax,q_stress,k_d_final,k_d_min_saved,k_d_max_saved,k_d_mean_saved,phi_c_final,clogging_amplification_final,L_final,L_ratio,Phi_final,Phi_ratio_total,R_post,tag",
    )
end

open(figure_csv, "w") do io
    println(
        io,
        "scenario,scenario_label,mode,Da,Re,BiC,stress_model,clogging_model,chi_H,phi_c_max,phi_c_ref_norm,sigma0_clean,sigma_eff_norm,kdmax,q_stress,k_d_final,phi_c_final,clogging_amplification_final,L_ratio,Phi_ratio_total,R_post",
    )
end

open(mapping_csv, "w") do io
    println(
        io,
        "scenario,scenario_label,Re,Ktilde,alpha_o,Sc_C,hb_over_lc,Sh_lam,Sh_turb,chi,Sh_eff,Bi,stress_model,clogging_model,chi_H,phi_c_max,phi_c_ref_norm,sigma0_clean,sigma_eff_ref_possible,sigma_eff_cap_possible,sigma_eff_norm",
    )
end

# ------------------------------------------------------------------------------
# Precompute Bi_C(Re), clean stresses and scenario-wise stress normalizations
# ------------------------------------------------------------------------------
maps = Dict{Tuple{Symbol,Float64},NamedTuple}()
stress_norms = Dict{Symbol,Float64}()

for pset in CFG.parameter_sets
    stress_norms[pset.name] = sigma_eff_norm_pset(pset)
end

for pset in CFG.parameter_sets, Re in pset.Re_values
    maps[(pset.name, Re)] = bi_from_re_pset(Re, pset)

    map = maps[(pset.name, Re)]
    clean = clean_stress_pset(Re, pset)
    sigma_eff_ref_possible = sigma_eff_ref_possible_pset(Re, pset)
    sigma_eff_cap_possible = sigma_eff_cap_possible_pset(Re, pset)

    open(mapping_csv, "a") do io
        append_csv_row(io, (
            String(pset.name),
            pset.label,
            Re,
            pset.Ktilde,
            pset.alpha_o,
            pset.Sc_C,
            pset.hb_over_lc,
            map.Sh_lam,
            map.Sh_turb,
            map.chi,
            map.Sh_eff,
            map.Bi,
            String(pset.stress_model),
            String(pset.clogging_model),
            pset.chi_H,
            pset.phi_c_max,
            pset.phi_c_ref_norm,
            clean.sigma0,
            sigma_eff_ref_possible,
            sigma_eff_cap_possible,
            stress_norms[pset.name],
        ))
    end
end

# ------------------------------------------------------------------------------
# Final sweep
# ------------------------------------------------------------------------------
@printf("\n=== Reynolds-coupled stress-based detachment sweep v2 ===\n")
@printf("mode = %s | antibiotic active\n", String(CFG.mode))
@printf("law: k_d(Re,L) = k_d,max * [sigma_eff(Re,L)/sigma_eff,ref]^q\n")
@printf("stress normalization = %s\n", String(CFG.stress_normalization))
@printf("chosen values: k_d,max = %.5f | q = %.3f\n", CFG.kdmax, CFG.q_stress)

for pset in CFG.parameter_sets, Da in pset.Da_values
    @printf("\n--- scenario = %s | Da = %.1f ---\n", String(pset.name), Da)

    local_rows = NamedTuple[]
    sigma_norm = stress_norms[pset.name]

    for Re in pset.Re_values
        map = maps[(pset.name, Re)]
        BiC = map.Bi
        kd_callback = make_kd_callback(Re, pset, sigma_norm)
        clean = clean_stress_pset(Re, pset)

        tag = @sprintf(
            "%s_Da_%0.1f_Re_%0.3e_stress_kdmax_%0.3f_q_%0.2f_chiH_%0.2f_phiref_%0.2f",
            String(pset.name),
            Da,
            Re,
            CFG.kdmax,
            CFG.q_stress,
            pset.chi_H,
            pset.phi_c_ref_norm,
        )

        _, _, summary = run_case_stress_detachment(;
            outdir = OUTDIR,
            tag = tag,

            L0tilde = CFG.L0tilde,
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
            Lmin = CFG.Lmin,

            # Nutrient-only Reynolds coupling, antibiotic active.
            C_top_bc = :robin,
            C_top_robin_Bi = BiC,
            C_top_robin_ref = 1.0,

            A_top_bc = :dirichlet,
            A_top = 0.0,
            A_step_time = CFG.tstart,
            A_step_value = CFG.A0,

            save_final_profiles = false,
            fixed_thickness = false,
            k_d_callback = kd_callback,
        )

        Phi0 = CFG.live_frac0 * CFG.L0tilde
        L_ratio = summary.L_final / CFG.L0tilde
        Phi_ratio_total = summary.Phi_final / max(Phi0, 1.0e-12)

        row = (
            scenario = String(pset.name),
            scenario_label = pset.label,
            mode = String(CFG.mode),
            Da = Da,
            Re = Re,
            Ktilde = pset.Ktilde,
            alpha_o = pset.alpha_o,
            Sc_C = pset.Sc_C,
            hb_over_lc = pset.hb_over_lc,
            BiC = BiC,
            Sh_eff = map.Sh_eff,
            Sh_lam = map.Sh_lam,
            Sh_turb = map.Sh_turb,
            chi = map.chi,
            stress_model = String(pset.stress_model),
            clogging_model = String(pset.clogging_model),
            chi_H = pset.chi_H,
            phi_c_max = pset.phi_c_max,
            phi_c_ref_norm = pset.phi_c_ref_norm,
            sigma0_clean = clean.sigma0,
            sigma_eff_norm = sigma_norm,
            kdmax = CFG.kdmax,
            q_stress = CFG.q_stress,
            k_d_final = summary.k_d_final,
            k_d_min_saved = safe_minimum(summary.kd_history),
            k_d_max_saved = safe_maximum(summary.kd_history),
            k_d_mean_saved = safe_mean(summary.kd_history),
            phi_c_final = summary.phi_c_final,
            clogging_amplification_final = summary.clogging_amplification_final,
            L_final = summary.L_final,
            L_ratio = L_ratio,
            Phi_final = summary.Phi_final,
            Phi_ratio_total = Phi_ratio_total,
            R_post = summary.R_final,
            tag = tag,
        )

        push!(local_rows, row)

        open(summary_csv, "a") do io
            append_csv_row(io, (
                row.scenario,
                row.scenario_label,
                row.mode,
                row.Da,
                row.Re,
                row.Ktilde,
                row.alpha_o,
                row.Sc_C,
                row.hb_over_lc,
                row.BiC,
                row.Sh_eff,
                row.Sh_lam,
                row.Sh_turb,
                row.chi,
                row.stress_model,
                row.clogging_model,
                row.chi_H,
                row.phi_c_max,
                row.phi_c_ref_norm,
                row.sigma0_clean,
                row.sigma_eff_norm,
                row.kdmax,
                row.q_stress,
                row.k_d_final,
                row.k_d_min_saved,
                row.k_d_max_saved,
                row.k_d_mean_saved,
                row.phi_c_final,
                row.clogging_amplification_final,
                row.L_final,
                row.L_ratio,
                row.Phi_final,
                row.Phi_ratio_total,
                row.R_post,
                row.tag,
            ))
        end

        open(figure_csv, "a") do io
            append_csv_row(io, (
                row.scenario,
                row.scenario_label,
                row.mode,
                row.Da,
                row.Re,
                row.BiC,
                row.stress_model,
                row.clogging_model,
                row.chi_H,
                row.phi_c_max,
                row.phi_c_ref_norm,
                row.sigma0_clean,
                row.sigma_eff_norm,
                row.kdmax,
                row.q_stress,
                row.k_d_final,
                row.phi_c_final,
                row.clogging_amplification_final,
                row.L_ratio,
                row.Phi_ratio_total,
                row.R_post,
            ))
        end

        @printf(
            "scenario = %s | Re = %.3e | BiC = %.4f | sigma0 = %.3e | sigma_norm = %.3e | kd_final = %.5f | phi_c = %.3f | amp = %.3f | L/L0 = %.4f | Phi/Phi0 = %.4f | R_post = %.4f\n",
            String(pset.name),
            Re,
            BiC,
            row.sigma0_clean,
            row.sigma_eff_norm,
            row.k_d_final,
            row.phi_c_final,
            row.clogging_amplification_final,
            row.L_ratio,
            row.Phi_ratio_total,
            row.R_post,
        )
    end

    Lvals = [r.L_ratio for r in local_rows]
    Pvals = [r.Phi_ratio_total for r in local_rows]
    Rvals = [r.R_post for r in local_rows]
    Kdvals = [r.k_d_final for r in local_rows]
    Phicvals = [r.phi_c_final for r in local_rows]
    Ampvals = [r.clogging_amplification_final for r in local_rows]

    Ltrend = classify_trend(Lvals)
    Ptrend = classify_trend(Pvals)
    Rtrend = classify_trend(Rvals)

    Re_vals = collect(pset.Re_values)
    Lpeak = peak_summary(Re_vals, Lvals)
    Ppeak = peak_summary(Re_vals, Pvals)
    Rpeak = peak_summary(Re_vals, Rvals)

    report_file = joinpath(
        OUTDIR,
        @sprintf("%s_Da_%0.1f_trend_report.txt", String(pset.name), Da),
    )

    open(report_file, "w") do io
        println(io, "Reynolds-coupled stress-based detachment sweep v2")
        println(io, "scenario = $(pset.name)")
        println(io, "scenario_label = $(pset.label)")
        println(io, "mode = $(CFG.mode)")
        println(io, "Da = $(Da)")
        println(io, "Ktilde = $(pset.Ktilde)")
        println(io, "alpha_o = $(pset.alpha_o)")
        println(io, "hb_over_lc = $(pset.hb_over_lc)")
        println(io, "stress_model = $(pset.stress_model)")
        println(io, "clogging_model = $(pset.clogging_model)")
        println(io, "chi_H = $(pset.chi_H)")
        println(io, "phi_c_max = $(pset.phi_c_max)")
        println(io, "phi_c_ref_norm = $(pset.phi_c_ref_norm)")
        println(io, "sigma_eff_norm = $(sigma_norm)")
        println(io, "k_d,max = $(CFG.kdmax)")
        println(io, "q = $(CFG.q_stress)")
        println(io, "law: k_d(Re,L) = k_d,max * [sigma_eff(Re,L)/sigma_eff,ref]^q")
        println(io, "")
        println(io, "Trend for L(T)/L(0):     " * Ltrend)
        println(io, "Trend for Phi(T)/Phi(0): " * Ptrend)
        println(io, "Trend for R_post:        " * Rtrend)
        println(io, "")
        println(io, "Peak summary for L(T)/L(0):")
        println(io, "  peak index = $(Lpeak.peak_idx)")
        println(io, "  peak Re    = $(Lpeak.peak_Re)")
        println(io, "  peak value = $(Lpeak.peak_val)")
        println(io, "")
        println(io, "Peak summary for Phi(T)/Phi(0):")
        println(io, "  peak index = $(Ppeak.peak_idx)")
        println(io, "  peak Re    = $(Ppeak.peak_Re)")
        println(io, "  peak value = $(Ppeak.peak_val)")
        println(io, "")
        println(io, "Peak summary for R_post:")
        println(io, "  peak index = $(Rpeak.peak_idx)")
        println(io, "  peak Re    = $(Rpeak.peak_Re)")
        println(io, "  peak value = $(Rpeak.peak_val)")
        println(io, "")
        println(io, "Re values = " * join(string.(pset.Re_values), ", "))
        println(io, "BiC values = " * join(string.([maps[(pset.name, Re)].Bi for Re in pset.Re_values]), ", "))
        println(io, "k_d_final values = " * join(string.(Kdvals), ", "))
        println(io, "phi_c_final values = " * join(string.(Phicvals), ", "))
        println(io, "clogging amplification values = " * join(string.(Ampvals), ", "))
        println(io, "L ratios = " * join(string.(Lvals), ", "))
        println(io, "Phi ratios = " * join(string.(Pvals), ", "))
        println(io, "R_post values = " * join(string.(Rvals), ", "))
    end
end

# ------------------------------------------------------------------------------
# Metadata note
# ------------------------------------------------------------------------------
meta_file = joinpath(OUTDIR, "README_re_stress_coupled_parameter_sets.txt")

open(meta_file, "w") do io
    println(io, "Reynolds-coupled stress-based detachment sweep v2")
    println(io, "------------------------------------------------")
    println(io, "Goal:")
    println(io, "  Generate Result 4 data after replacing the previous Bi_C-based")
    println(io, "  detachment law by a stress-based effective closure.")
    println(io, "")
    println(io, "Closure:")
    println(io, "  k_d(Re,L) = k_d,max * [sigma_eff(Re,L)/sigma_eff,ref]^q")
    println(io, "  k_d,max = $(CFG.kdmax)")
    println(io, "  q       = $(CFG.q_stress)")
    println(io, "  stress_normalization = $(CFG.stress_normalization)")
    println(io, "")
    println(io, "Important correction:")
    println(io, "  phi_c_max is only a numerical cap.")
    println(io, "  The stress normalization uses phi_c_ref_norm.")
    println(io, "  This avoids over-normalising sigma_eff by the near-singular")
    println(io, "  amplification associated with phi_c_max.")
    println(io, "")
    println(io, "Physical interpretation:")
    println(io, "  Bi_C(Re) controls external nutrient transfer.")
    println(io, "  sigma_eff(Re,L) controls mechanical removal.")
    println(io, "  In the microfluidic-effective case, sigma_eff includes a clogging")
    println(io, "  amplification factor (1 - phi_c)^(-2).")
    println(io, "")
    println(io, "Chosen setup:")
    println(io, "  mode = nutrient_only")
    println(io, "  antibiotic remains active")
    println(io, "  A0 = $(CFG.A0)")
    println(io, "")
    println(io, "Common settings:")
    println(io, "  L0tilde  = $(CFG.L0tilde)")
    println(io, "  live0    = $(CFG.live_frac0)")
    println(io, "  Tfinal   = $(CFG.Tfinal)")
    println(io, "  tstart   = $(CFG.tstart)")
    println(io, "  Lambda   = $(CFG.Lambda)")
    println(io, "  beta_Ab  = $(CFG.beta_Ab)")
    println(io, "  Da_a     = $(CFG.Da_a)")
    println(io, "  PiE0     = $(CFG.PiE0)")
    println(io, "  p        = $(CFG.kill_power)")
    println(io, "  gamma    = $(CFG.gamma_kill)")
    println(io, "")
    println(io, "Parameter sets:")

    for pset in CFG.parameter_sets
        println(io, "  $(pset.name): $(pset.label)")
        println(io, "    Ktilde = $(pset.Ktilde)")
        println(io, "    alpha_o = $(pset.alpha_o)")
        println(io, "    Sc_C = $(pset.Sc_C)")
        println(io, "    hb_over_lc = $(pset.hb_over_lc)")
        println(io, "    Da values = " * join(string.(pset.Da_values), ", "))
        println(io, "    Re values = " * join(string.(pset.Re_values), ", "))
        println(io, "    stress_model = $(pset.stress_model)")
        println(io, "    clogging_model = $(pset.clogging_model)")
        println(io, "    chi_H = $(pset.chi_H)")
        println(io, "    phi_c_max = $(pset.phi_c_max)")
        println(io, "    phi_c_ref_norm = $(pset.phi_c_ref_norm)")
        println(io, "    sigma_eff_norm = $(stress_norms[pset.name])")
    end

    println(io, "")
    println(io, "Main files:")
    println(io, "  re_stress_coupled_parameter_sets_summary.csv")
    println(io, "  re_stress_coupled_parameter_sets_figure_ready.csv")
    println(io, "  re_stress_coupled_parameter_sets_mapping.csv")
    println(io, "  one trend report per scenario and Da")
end

println("\nDone.")
println("Output folder:    ", OUTDIR)
println("Summary CSV:      ", summary_csv)
println("Figure-ready CSV: ", figure_csv)
println("Mapping CSV:      ", mapping_csv)
println("Metadata note:    ", meta_file)