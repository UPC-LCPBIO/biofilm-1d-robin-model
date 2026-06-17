module Grid1D

export UniformGrid1D, make_grid, trapz_uniform

struct UniformGrid1D
    H::Float64
    N::Int
    z::Vector{Float64}
    dz::Float64
end

function make_grid(H::Real, N::Integer)
    @assert N ≥ 3
    z = collect(range(0.0, Float64(H), length=N))
    dz = z[2] - z[1]
    return UniformGrid1D(Float64(H), Int(N), z, dz)
end


function trapz_uniform(f::AbstractVector{<:Real}, dz::Real)
    n = length(f)
    @assert n ≥ 2
    s = 0.5*(f[1] + f[end])
    for i in 2:n-1
        s += f[i]
    end
    return Float64(dz) * Float64(s)
end

end
