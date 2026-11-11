function reference_diffusion_profile(z, t, L; nterms=100)
    total = 0.0
    for k in 0:nterms-1
        j = 2k + 1
        lambda = j * pi / (2L)
        total += 4 * (-1)^k / (j * pi) * cos(lambda * z) * exp(-lambda^2 * t)
    end
    return 1.0 - total
end

@testset "Diffusion grid refinement" begin
    errors = Float64[]
    for n in (17, 33, 65)
        z, _, s = run_case(;
            N=n, Tfinal=0.1, dt=5e-5, Lambda=0.0, Da=0.0,
            C_init=0.0, C_top_bc=:dirichlet, C_top=1.0,
            fixed_thickness=true, write_output=false)
        exact = [reference_diffusion_profile(zi, 0.1, 0.8) for zi in z]
        push!(errors, sqrt(sum(abs2, s.C_final .- exact) / n))
    end
    @test errors[2] < errors[1]
    @test errors[3] < errors[2]
end
