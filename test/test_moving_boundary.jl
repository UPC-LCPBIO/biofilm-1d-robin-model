@testset "Moving mesh geometric conservation" begin
    # Uniform, fully live biofilm with constant nutrient remains phi == 1
    # while growth increases thickness. This tests the ALE Jacobian and flux.
    _, L, s = run_case(;
        N=41, Tfinal=0.2, dt=0.001, Lambda=0.2,
        Da=0.0, C_init=1.0, C_top_bc=:dirichlet, C_top=1.0,
        live_frac0=1.0, A_step_time=1.0, A_top_bc=:dirichlet,
        write_output=false)
    @test L > 0.8
    @test maximum(abs.(s.phi_final .- 1.0)) < 1e-10
    @test maximum(abs.(s.C_final .- 1.0)) < 1e-10
    @test abs(s.conservation_residual) < 1e-10

    # Pure surface detachment removes volume, not the surviving volume fraction.
    _, Ldet, detached = run_case(;
        N=41, Tfinal=0.2, dt=0.001, Lambda=0.0,
        Da=0.0, C_init=0.0, k_d=0.05, write_output=false)
    @test Ldet < 0.8
    @test maximum(abs.(detached.phi_final .- 0.6)) < 1e-10
    @test detached.detachment_integrated > 0
    @test abs(detached.conservation_residual) < 1e-10
end
