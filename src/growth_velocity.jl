"""Integrate v_z = Lambda * e(C) * phi, with v(0) = 0."""
function compute_growth_velocity!(v, e, phi, Lambda::Float64, L::Float64, dxi::Float64)
    v[1] = 0.0
    @inbounds for i in 2:length(v)
        v[i] = v[i - 1] + 0.5 * Lambda * L * dxi *
            (e[i - 1] * phi[i - 1] + e[i] * phi[i])
    end
    return v
end
