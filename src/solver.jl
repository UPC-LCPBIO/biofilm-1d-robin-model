function write_timeseries(path::String, t, L, Phi, vtop, B, Cavg, Aavg, kd)
    open(path, "w") do io
        println(io, "t,L,Phi,v_top,B,Cavg,Aavg,k_d")
        for i in eachindex(t)
            @printf(io, "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n",
                    t[i], L[i], Phi[i], vtop[i], B[i], Cavg[i], Aavg[i], kd[i])
        end
    end
    return path
end

function write_profiles(path::String, z, C, A, phi, v, e)
    open(path, "w") do io
        println(io, "z,C,A,phi_l,v,e")
        for i in eachindex(z)
            @printf(io, "%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n",
                    z[i], C[i], A[i], phi[i], v[i], e[i])
        end
    end
    return path
end

"""
    run_case(; kwargs...)

Integrate the one-dimensional biofilm equations on a fixed reference grid.
Biomass is updated by conservative ALE finite-volume fluxes. Solute diffusion
is implicit, while reaction and geometric advection are explicit.

`k_d_law(t, L)` returns the coefficient in `k_d * L^2`.
`fixed_thickness=true` retains the original no-expansion, no-advection model.
The output tuple and standard CSV files retain the paper scripts' interfaces.
"""
function run_case(;
    L0tilde::Float64=L0tilde,
    N::Int=N,
    Tfinal::Float64=Tfinal,
    dt::Float64=dt,
    save_every::Int=save_every,
    Lambda::Float64=Lambda,
    Da::Float64=Da,
    Da_a::Float64=Da_a,
    beta_Ab::Float64=beta_Ab,
    Ktilde::Float64=Ktilde,
    PiE0::Float64=PiE0,
    kill_power::Float64=kill_power,
    gamma_kill::Float64=gamma_kill,
    live_frac0::Float64=live_frac0,
    k_d::Float64=k_d,
    k_d_law::Union{Nothing,Function}=nothing,
    Lmin::Float64=Lmin,
    C_init::Float64=C_init,
    A_init::Float64=A_init,
    C_top_bc::Symbol=C_top_bc,
    C_top::Float64=C_top,
    C_top_robin_Bi::Float64=C_top_robin_Bi,
    C_top_robin_ref::Float64=C_top_robin_ref,
    A_top_bc::Symbol=A_top_bc,
    A_top::Float64=A_top,
    A_top_robin_Bi::Float64=A_top_robin_Bi,
    A_step_time::Float64=A_step_time,
    A_step_value::Float64=A_step_value,
    outdir::Union{Nothing,String}="output",
    tag::String="case",
    save_final_profiles::Bool=true,
    write_output::Bool=true,
    fixed_thickness::Bool=false,
)
    C_top_bc in (:dirichlet, :robin) || throw(ArgumentError("Invalid nutrient boundary condition"))
    A_top_bc in (:dirichlet, :robin) || throw(ArgumentError("Invalid antibiotic boundary condition"))
    0 < live_frac0 <= 1 || throw(ArgumentError("Invalid initial live fraction"))
    N >= 3 || throw(ArgumentError("At least three grid nodes are required"))
    dt > 0 && Tfinal > 0 && L0tilde > 0 || throw(ArgumentError("Invalid time or length"))
    save_every >= 1 || throw(ArgumentError("save_every must be positive"))
    beta_Ab > 0 && Ktilde > 0 && Lmin > 0 || throw(ArgumentError("Invalid diffusivity or length"))
    C_top_robin_Bi >= 0 && A_top_robin_Bi >= 0 || throw(ArgumentError("Biot numbers must be nonnegative"))
    k_d >= 0 || throw(ArgumentError("Detachment coefficient must be nonnegative"))
    gamma_kill > 0 && kill_power >= 0 || throw(ArgumentError("Invalid killing exponents"))

    grid = make_reference_grid(N)
    L = L0tilde
    C = fill(C_init, N)
    A = fill(A_init, N)
    phi = fill(live_frac0, N)
    e = zeros(N)
    v = zeros(N)
    RC = zeros(N)
    RA = zeros(N)
    work = TransportWorkspace(N)

    t_hist = Float64[]
    L_hist = Float64[]
    Phi_hist = Float64[]
    vtop_hist = Float64[]
    B_hist = Float64[]
    Cavg_hist = Float64[]
    Aavg_hist = Float64[]
    kd_hist = Float64[]

    Phi_initial = L * weighted_sum(phi, grid.weights)
    Phi_tstart = NaN
    G_integrated = 0.0
    K_integrated = 0.0
    D_integrated = 0.0
    clipped_biomass = 0.0
    clipped_solute = 0.0
    max_live_fraction = maximum(phi)
    min_live_fraction = minimum(phi)
    t = 0.0
    istep = 0

    if write_output && outdir !== nothing
        mkpath(outdir)
    end

    while t < Tfinal - 1e-14
        istep += 1
        step = min(dt, Tfinal - t)
        tnew = t + step

        monod!(e, C, Ktilde)
        if fixed_thickness
            fill!(v, 0.0)
        else
            compute_growth_velocity!(v, e, phi, Lambda, L, grid.dxi)
        end

        kd_now = fixed_thickness ? 0.0 :
            (k_d_law === nothing ? k_d : Float64(k_d_law(t, L)))
        isfinite(kd_now) && kd_now >= 0 || throw(DomainError(kd_now, "Invalid detachment coefficient"))
        detachment_velocity = kd_now * L^2
        Ldot = fixed_thickness ? 0.0 : v[end] - detachment_velocity
        Lnew = fixed_thickness ? L : max(L + step * Ldot, Lmin)
        Ldot = (Lnew - L) / step

        @inbounds for i in 1:N
            RC[i] = -Da * e[i] * phi[i]
            RA[i] = -Da_a * A[i] * phi[i]
        end

        nutrient_ref = C_top_bc === :robin ? C_top_robin_ref : C_top
        clipped_solute += step_solute_ale!(C, grid, work, L, Lnew, Ldot,
            step, 1.0, RC, C_top_bc, nutrient_ref, C_top_robin_Bi)

        antibiotic_ref = tnew >= A_step_time ? A_step_value : A_top
        clipped_solute += step_solute_ale!(A, grid, work, L, Lnew, Ldot,
            step, beta_Ab, RA, A_top_bc, antibiotic_ref, A_top_robin_Bi)

        growth, killing, detached, clipped = step_biomass_ale!(phi, e, A, v,
            grid, work, L, Lnew, Ldot, step, Lambda, PiE0, kill_power, gamma_kill)
        G_integrated += step * growth
        K_integrated += step * killing
        D_integrated += step * detached
        clipped_biomass += clipped

        L = Lnew
        t = tnew
        max_live_fraction = max(max_live_fraction, maximum(phi))
        min_live_fraction = min(min_live_fraction, minimum(phi))
        if isnan(Phi_tstart) && t >= A_step_time
            Phi_tstart = L * weighted_sum(phi, grid.weights)
        end

        if istep % save_every == 0 || t >= Tfinal - 1e-14
            Phi = L * weighted_sum(phi, grid.weights)
            Cavg = weighted_sum(C, grid.weights)
            Aavg = weighted_sum(A, grid.weights)
            B = isnan(Phi_tstart) ? NaN : Phi / max(Phi_tstart, 1e-12)
            kd_save = fixed_thickness ? 0.0 :
                (k_d_law === nothing ? k_d : Float64(k_d_law(t, L)))

            push!(t_hist, t)
            push!(L_hist, L)
            push!(Phi_hist, Phi)
            push!(vtop_hist, v[end])
            push!(B_hist, B)
            push!(Cavg_hist, Cavg)
            push!(Aavg_hist, Aavg)
            push!(kd_hist, kd_save)
        end
    end

    monod!(e, C, Ktilde)
    if fixed_thickness
        fill!(v, 0.0)
    else
        compute_growth_velocity!(v, e, phi, Lambda, L, grid.dxi)
    end
    z = grid.xi .* L
    Phi_final = L * weighted_sum(phi, grid.weights)
    R_final = isnan(Phi_tstart) ? NaN : Phi_final / max(Phi_tstart, 1e-12)
    conservation_residual = Phi_final - Phi_initial -
        (G_integrated - K_integrated - D_integrated + clipped_biomass)

    summary = (
        L_final=L,
        Phi_final=Phi_final,
        Phi_tstart=Phi_tstart,
        R_final=R_final,
        G_final=L / L0tilde,
        t_history=t_hist,
        L_history=L_hist,
        Phi_history=Phi_hist,
        vtop_history=vtop_hist,
        B_history=B_hist,
        Cavg_history=Cavg_hist,
        Aavg_history=Aavg_hist,
        kd_history=kd_hist,
        z=z,
        C_final=copy(C),
        A_final=copy(A),
        phi_final=copy(phi),
        v_final=copy(v),
        e_final=copy(e),
        Phi_initial=Phi_initial,
        growth_integrated=G_integrated,
        killing_integrated=K_integrated,
        detachment_integrated=D_integrated,
        conservation_residual=conservation_residual,
        clipped_biomass=clipped_biomass,
        clipped_solute=clipped_solute,
        max_live_fraction=max_live_fraction,
        min_live_fraction=min_live_fraction,
    )

    if write_output && outdir !== nothing
        write_timeseries(joinpath(outdir, "$(tag)_timeseries.csv"),
            t_hist, L_hist, Phi_hist, vtop_hist, B_hist, Cavg_hist, Aavg_hist, kd_hist)
        if save_final_profiles
            write_profiles(joinpath(outdir, "$(tag)_final_profiles.csv"), z, C, A, phi, v, e)
        end
    end
    return z, L, summary
end
