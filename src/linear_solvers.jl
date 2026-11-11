"""In-place Thomas solver; arrays use a[i] for the lower diagonal."""
function thomas_solve!(x, a, b, c, rhs, cprime, dprime)
    n = length(b)
    @assert length(a) == n && length(c) == n && length(rhs) == n
    @assert length(x) == n && length(cprime) == n && length(dprime) == n
    @inbounds begin
        cprime[1] = c[1] / b[1]
        dprime[1] = rhs[1] / b[1]
        for i in 2:n
            pivot = b[i] - a[i] * cprime[i - 1]
            cprime[i] = c[i] / pivot
            dprime[i] = (rhs[i] - a[i] * dprime[i - 1]) / pivot
        end
        x[n] = dprime[n]
        for i in (n - 1):-1:1
            x[i] = dprime[i] - cprime[i] * x[i + 1]
        end
    end
    return x
end
