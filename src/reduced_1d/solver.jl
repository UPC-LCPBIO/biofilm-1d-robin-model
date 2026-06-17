include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "grid.jl"))
include(joinpath(@__DIR__, "linear_solvers.jl"))
include(joinpath(@__DIR__, "transport.jl"))
include(joinpath(@__DIR__, "growth_velocity.jl"))

using .Grid1D
using .Transport1D
using .ConstraintNu1D
using DelimitedFiles
using Printf

function linear_remap(z_old::Vector{Float64}, f_old::Vector{Float64}, z_new::Vector{Float64})
    out = similar(z_new)
    j = 1
    nold = length(z_old)
    @inbounds for i in eachindex(z_new)
        x = z_new[i]
        if x <= z_old[1]
            out[i] = f_old[1]
        elseif x >= z_old[end]
            out[i] = f_old[end]
        else
            while j < nold - 1 && z_old[j + 1] < x
                j += 1
            end
            x0, x1 = z_old[j], z_old[j + 1]
            y0, y1 = f_old[j], f_old[j + 1]
            θ = (x - x0) / (x1 - x0)
            out[i] = (1 - θ) * y0 + θ * y1
        end
    end
    return out
end

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
    @assert C_top_bc in (:dirichlet, :robin)
    @assert A_top_bc in (:dirichlet, :robin)
    @assert 0.0 < live_frac0 <= 1.0
    @assert N ≥ 3
    @assert dt > 0.0 && Tfinal > 0.0
    @assert save_every ≥ 1

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

    Phi_tstart = NaN
    t = 0.0
    istep = 0

    if write_output && outdir !== nothing
        mkpath(outdir)
    end

    while t < Tfinal - 1.0e-15
        istep += 1
        Δt = min(dt, Tfinal - t)
        t += Δt

        Transport1D.monod!(e, C, Ktilde)

        if fixed_thickness
            fill!(v, 0.0)
        else
            ConstraintNu1D.compute_v!(v, e, phi, Lambda, dz)
        end

        @inbounds for i in 1:N
            R_C[i] = -Da * e[i] * phi[i]
            R_A[i] = -Da_a * A[i] * phi[i]
            R_phi[i] = Lambda * e[i] * phi[i] - Lambda * PiE0 * (e[i]^kill_power) * A[i] * (phi[i]^gamma_kill)
        end

        if C_top_bc === :robin
            Transport1D.step_scalar_imex!(C, v_zero, D_C, Δt, dz;
                reaction_rhs=R_C,
                top_dirichlet=false,
                top_robin=true,
                top_robin_coeff=C_top_robin_Bi,
                top_robin_ref=C_top_robin_ref,
            )
        else
            Transport1D.step_scalar_imex!(C, v_zero, D_C, Δt, dz;
                reaction_rhs=R_C,
                top_dirichlet=true,
                top_value=C_top,
            )
        end

        A_bc_now = (t >= A_step_time) ? A_step_value : A_top
        if A_top_bc === :robin
            Transport1D.step_scalar_imex!(A, v_zero, D_A, Δt, dz;
                reaction_rhs=R_A,
                top_dirichlet=false,
                top_robin=true,
                top_robin_coeff=A_top_robin_Bi,
                top_robin_ref=A_bc_now,
            )
        else
            Transport1D.step_scalar_imex!(A, v_zero, D_A, Δt, dz;
                reaction_rhs=R_A,
                top_dirichlet=true,
                top_value=A_bc_now,
            )
        end

        Transport1D.step_scalar_imex!(phi, v, D_zero, Δt, dz;
            reaction_rhs=R_phi,
            top_dirichlet=false,
            advect_top_outflow=true,
        )

        @inbounds for i in 1:N
            C[i] = max(C[i], 0.0)
            A[i] = max(A[i], 0.0)
            phi[i] = max(phi[i], 0.0)
        end

        kd_now = k_d_law === nothing ? k_d : Float64(k_d_law(t, L))

        if !fixed_thickness
            Ldot = v[end] - kd_now * L^2
            Lnew = max(L + Δt * Ldot, Lmin)

            if abs(Lnew - L) > 1.0e-14
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

        if istep % save_every == 0 || t >= Tfinal - 1.0e-15
            Phi = Grid1D.trapz_uniform(phi, dz)
            Cavg = Grid1D.trapz_uniform(C, dz) / L
            Aavg = Grid1D.trapz_uniform(A, dz) / L
            B = isnan(Phi_tstart) ? NaN : (Phi / max(Phi_tstart, 1.0e-12))

            push!(t_hist, t)
            push!(L_hist, L)
            push!(Phi_hist, Phi)
            push!(vtop_hist, v[end])
            push!(B_hist, B)
            push!(Cavg_hist, Cavg)
            push!(Aavg_hist, Aavg)
            push!(kd_hist, kd_now)
        end
    end

    Phi_final = Grid1D.trapz_uniform(phi, dz)
    R_final = isnan(Phi_tstart) ? NaN : (Phi_final / max(Phi_tstart, 1.0e-12))

    summary = (
        L_final = L,
        Phi_final = Phi_final,
        Phi_tstart = Phi_tstart,
        R_final = R_final,
        G_final = L / L0tilde,
        t_history = t_hist,
        L_history = L_hist,
        Phi_history = Phi_hist,
        vtop_history = vtop_hist,
        B_history = B_hist,
        Cavg_history = Cavg_hist,
        Aavg_history = Aavg_hist,
        kd_history = kd_hist,
        z = z,
        C_final = copy(C),
        A_final = copy(A),
        phi_final = copy(phi),
        v_final = copy(v),
        e_final = copy(e),
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
