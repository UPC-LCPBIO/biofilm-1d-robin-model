@testset "Public result interface and detachment callback" begin
    callback = (t, L) -> 0.01 * (1.0 + 0.1 * t)
    mktempdir() do dir
        z, L, s = run_case(;
            N=33, dt=0.001, Tfinal=0.04, save_every=10,
            Lambda=0.2, Da=6.0, k_d_law=callback,
            outdir=dir, tag="smoke", write_output=true)
        @test isfile(joinpath(dir, "smoke_timeseries.csv"))
        @test isfile(joinpath(dir, "smoke_final_profiles.csv"))
        @test L == s.L_final
        @test z == s.z
        @test length(s.t_history) == 4
        @test s.kd_history[end] ≈ callback(0.04, L) atol=1e-12
        @test all(isfinite, s.Phi_history)
        @test isfinite(s.G_final) && isnan(s.R_final)
        @test abs(s.conservation_residual) < 1e-10
    end

    # Fixed-thickness growth reduces to an explicit Euler scalar ODE when
    # nutrient is uniform and there is no diffusion-reaction coupling.
    _, _, growth = run_case(;
        N=33, Tfinal=0.2, dt=0.001, Lambda=0.2, Da=0.0,
        C_init=1.0, C_top_bc=:dirichlet, C_top=1.0,
        fixed_thickness=true, write_output=false)
    econstant = 1.0 / (1.0 + 2.857142857142857)
    expected_phi = 0.6 * (1.0 + 0.001 * 0.2 * econstant)^200
    @test maximum(abs.(growth.phi_final .- expected_phi)) < 1e-10

    _, Lfix, fixed = run_case(;
        N=33, dt=0.001, Tfinal=0.2, Lambda=0.0,
        k_d=0.05, write_output=false, fixed_thickness=true)
    @test Lfix == 0.8
    @test fixed.detachment_integrated == 0.0
    @test fixed.growth_integrated == 0.0
end
