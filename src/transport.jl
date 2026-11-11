"""Work arrays reused by the conservative ALE transport steps."""
struct TransportWorkspace
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    rhs::Vector{Float64}
    cp::Vector{Float64}
    dp::Vector{Float64}
    x::Vector{Float64}
    flux::Vector{Float64}
    amount::Vector{Float64}
end

function TransportWorkspace(N::Int)
    return TransportWorkspace(
        zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N),
        zeros(N), zeros(N + 1), zeros(N))
end

function monod!(out, C, K::Float64)
    @inbounds for i in eachindex(C)
        out[i] = C[i] / (C[i] + K)
    end
    return out
end

"""
Advance a solute on the reference domain with explicit geometric flux and
reaction, implicit physical diffusion, and a top Robin or Dirichlet boundary.
The physical solute velocity is zero in the original model.
"""
function step_solute_ale!(S, grid::ReferenceGrid1D, work::TransportWorkspace,
                          L::Float64, Lnew::Float64, Ldot::Float64,
                          dt::Float64, diffusivity::Float64, reaction,
                          mode::Symbol, top_value::Float64, Bi::Float64)
    N = grid.N
    xi, weights, dxi = grid.xi, grid.weights, grid.dxi
    flux = work.flux
    flux[1] = 0.0
    @inbounds for k in 2:N
        xface = 0.5 * (xi[k - 1] + xi[k])
        u = -xface * Ldot
        flux[k] = u * (u >= 0 ? S[k - 1] : S[k])
    end
    flux[N + 1] = -Ldot * S[N]

    kappa = diffusivity / (Lnew * dxi)
    @inbounds for i in 1:N
        work.rhs[i] = weights[i] * L * (S[i] + dt * reaction[i]) -
                      dt * (flux[i + 1] - flux[i])
        work.a[i] = 0.0
        work.c[i] = 0.0
        work.b[i] = weights[i] * Lnew
        if i > 1
            work.a[i] = -dt * kappa
            work.b[i] += dt * kappa
        end
        if i < N
            work.c[i] = -dt * kappa
            work.b[i] += dt * kappa
        end
    end
    if mode === :robin
        work.b[N] += dt * Bi
        work.rhs[N] += dt * Bi * top_value
    elseif mode === :dirichlet
        work.a[N] = 0.0
        work.b[N] = 1.0
        work.c[N] = 0.0
        work.rhs[N] = top_value
    else
        throw(ArgumentError("Unknown top boundary condition: $mode"))
    end
    thomas_solve!(work.x, work.a, work.b, work.c, work.rhs, work.cp, work.dp)
    removed = 0.0
    @inbounds for i in 1:N
        if work.x[i] < 0.0
            removed += -work.x[i] * weights[i] * Lnew
            work.x[i] = 0.0
        end
        S[i] = work.x[i]
    end
    return removed
end

"""
Advance live biomass as cell amounts L*w[i]*phi[i]. The upper boundary flux
is the detachment velocity relative to the expanding reference mesh.
"""
function step_biomass_ale!(phi, e, A, v, grid::ReferenceGrid1D,
                           work::TransportWorkspace, L::Float64, Lnew::Float64,
                           Ldot::Float64, dt::Float64, Lambda::Float64,
                           PiE0::Float64, kill_power::Float64, gamma_kill::Float64)
    N = grid.N
    xi, w = grid.xi, grid.weights
    F = work.flux
    F[1] = 0.0
    @inbounds for k in 2:N
        xface = 0.5 * (xi[k - 1] + xi[k])
        U = 0.5 * (v[k - 1] + v[k]) - xface * Ldot
        F[k] = U * (U >= 0.0 ? phi[k - 1] : phi[k])
    end
    F[N + 1] = (v[N] - Ldot) * phi[N]
    growth = 0.0
    killing = 0.0
    detachment = F[N + 1]
    @inbounds for i in 1:N
        gi = Lambda * e[i] * phi[i]
        ki = Lambda * PiE0 * e[i]^kill_power * A[i] * max(phi[i], 0.0)^gamma_kill
        growth += L * w[i] * gi
        killing += L * w[i] * ki
        work.amount[i] = L * w[i] * phi[i] - dt * (F[i + 1] - F[i]) +
                         dt * L * w[i] * (gi - ki)
    end
    clipped = 0.0
    @inbounds for i in 1:N
        if work.amount[i] < 0.0
            clipped -= work.amount[i]
            work.amount[i] = 0.0
        end
        phi[i] = work.amount[i] / (Lnew * w[i])
    end
    return growth, killing, detachment, clipped
end
