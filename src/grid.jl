"""Uniform reference grid with node-centered control volumes."""
struct ReferenceGrid1D
    N::Int
    xi::Vector{Float64}
    weights::Vector{Float64}
    dxi::Float64
end

function make_reference_grid(N::Integer)
    N >= 3 || throw(ArgumentError("At least three grid nodes are required"))
    dxi = 1.0 / (N - 1)
    xi = collect(range(0.0, 1.0; length=N))
    weights = fill(dxi, N)
    weights[1] = weights[end] = dxi / 2
    return ReferenceGrid1D(Int(N), xi, weights, dxi)
end

function weighted_sum(values::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    length(values) == length(weights) || throw(DimensionMismatch("Quadrature weights"))
    acc = 0.0
    @inbounds for i in eachindex(values, weights)
        acc += weights[i] * values[i]
    end
    return acc
end
