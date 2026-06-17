module Transport1D

using ..Solvers1D: build_diffusion_tridiag, thomas_solve!, apply_dirichlet_top!

export monod!, upwind_flux_faces!, step_scalar_imex!

function monod!(out::Vector{Float64}, C::Vector{Float64}, K::Float64)
    N = length(C)
    @assert length(out) == N
    @inbounds for i in 1:N
        out[i] = C[i] / (C[i] + K)
    end
    return out
end

function upwind_flux_faces!(Ffaces::Vector{Float64}, v::Vector{Float64}, s::Vector{Float64})
    N = length(s)
    @assert length(v) == N
    @assert length(Ffaces) == N - 1
    @inbounds for i in 1:N-1
        vface = 0.5 * (v[i] + v[i+1])
        Ffaces[i] = (vface >= 0.0) ? (vface * s[i]) : (vface * s[i+1])
    end
    return Ffaces
end

"""
IMEX step for
    s_t = -(v s)_z + (D s_z)_z + reaction_rhs
on a uniform 1D grid with default no-flux at z=0.
At the top boundary, either Dirichlet or Robin can be imposed.
"""
function step_scalar_imex!(
    s::Vector{Float64},
    v::Vector{Float64},
    D::Vector{Float64},
    dt::Float64,
    dz::Float64;
    reaction_rhs::Vector{Float64},
    top_dirichlet::Bool=true,
    top_value::Float64=1.0,
    top_robin::Bool=false,
    top_robin_coeff::Float64=0.0,
    top_robin_ref::Float64=1.0,
    advect_top_outflow::Bool=false,
)
    N = length(s)
    @assert N >= 2
    @assert length(v) == N && length(D) == N && length(reaction_rhs) == N
    @assert !(top_dirichlet && top_robin) "Use either top Dirichlet or top Robin, not both."

    Ffaces = zeros(Float64, N-1)
    upwind_flux_faces!(Ffaces, v, s)

    adv = zeros(Float64, N)
    invdz = 1.0 / dz
    Fbot = 0.0

    Ftop = 0.0
    if advect_top_outflow
        vN = v[N]
        Ftop = (vN >= 0.0) ? (vN * s[N]) : 0.0
    end

    adv[1] = -(Ffaces[1] - Fbot) * invdz
    @inbounds for i in 2:N-1
        adv[i] = -(Ffaces[i] - Ffaces[i-1]) * invdz
    end
    adv[N] = -(Ftop - Ffaces[N-1]) * invdz

    rhs = similar(s)
    @inbounds for i in 1:N
        rhs[i] = s[i] + dt * (adv[i] + reaction_rhs[i])
    end

    a, b, c = build_diffusion_tridiag(D, dt, dz)

    if top_robin
        Dface = 0.5 * (D[N-1] + D[N])
        a[N] = -Dface / dz
        b[N] =  Dface / dz + top_robin_coeff
        c[N] = 0.0
        rhs[N] = top_robin_coeff * top_robin_ref
    elseif top_dirichlet
        apply_dirichlet_top!(a, b, c, rhs, top_value)
    end

    thomas_solve!(a, b, c, rhs)
    copyto!(s, rhs)
    return s
end

end
