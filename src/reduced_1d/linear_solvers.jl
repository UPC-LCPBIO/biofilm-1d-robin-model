module Solvers1D

export thomas_solve!, build_diffusion_tridiag, apply_dirichlet_top!, apply_dirichlet_bottom!

function thomas_solve!(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, rhs::Vector{Float64})
    n = length(b)
    @assert length(a)==n && length(c)==n && length(rhs)==n

    for i in 2:n
        w = a[i]/b[i-1]
        b[i] -= w*c[i-1]
        rhs[i] -= w*rhs[i-1]
    end

    rhs[n] /= b[n]
    for i in (n-1):-1:1
        rhs[i] = (rhs[i] - c[i]*rhs[i+1]) / b[i]
    end
    return rhs
end


function build_diffusion_tridiag(D::Vector{Float64}, dt::Float64, dz::Float64)
    N = length(D)
    a = zeros(Float64, N)
    b = ones(Float64, N)
    c = zeros(Float64, N)

    invdz2 = 1.0/(dz*dz)

    for i in 2:N-1
        Dm = 0.5*(D[i-1] + D[i])
        Dp = 0.5*(D[i] + D[i+1])
        a[i] = -dt * Dm * invdz2
        c[i] = -dt * Dp * invdz2
        b[i] = 1.0 - (a[i] + c[i])
    end

    Dp = 0.5*(D[1] + D[2])
    c[1] = -dt * Dp * invdz2
    b[1] = 1.0 - c[1]
    a[1] = 0.0

    Dm = 0.5*(D[N-1] + D[N])
    a[N] = -dt * Dm * invdz2
    b[N] = 1.0 - a[N]
    c[N] = 0.0

    return a, b, c
end

function apply_dirichlet_top!(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, rhs::Vector{Float64}, value::Float64)
    N = length(b)
    a[N] = 0.0
    b[N] = 1.0
    c[N] = 0.0
    rhs[N] = value
    return nothing
end

function apply_dirichlet_bottom!(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, rhs::Vector{Float64}, value::Float64)
    a[1] = 0.0
    b[1] = 1.0
    c[1] = 0.0
    rhs[1] = value
    return nothing
end

end
