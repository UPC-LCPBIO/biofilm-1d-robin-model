module ConstraintNu1D

export compute_v!

"""
    compute_v!(v, e, phi, Lambda, dz)

Computes the growth-induced velocity field from
    ∂_z v = Λ e φ,   v(0)=0
using trapezoidal integration on a uniform grid.
"""
function compute_v!(v::Vector{Float64}, e::Vector{Float64}, phi::Vector{Float64}, Lambda::Float64, dz::Float64)
    N = length(v)
    @assert length(e) == N && length(phi) == N

    v[1] = 0.0
    acc = 0.0
    @inbounds for i in 2:N
        rhs_im1 = Lambda * e[i-1] * phi[i-1]
        rhs_i   = Lambda * e[i]   * phi[i]
        acc += 0.5 * (rhs_im1 + rhs_i) * dz
        v[i] = acc
    end
    return v
end

end
