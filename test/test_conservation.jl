@testset "Conservative biomass budget" begin
    _, _, s = run_case(;
        N=65, dt=2e-4, Tfinal=6.1, save_every=1000,
        Da=18.0, C_top_bc=:robin, C_top_robin_Bi=5.0,
        A_top_bc=:dirichlet, A_step_time=0.5, A_step_value=0.05,
        k_d=0.01, write_output=false)

    @test isfinite(s.Phi_final)
    @test abs(s.conservation_residual) < 1e-9
    @test abs(s.Phi_final - s.Phi_initial -
        (s.growth_integrated - s.killing_integrated -
         s.detachment_integrated + s.clipped_biomass)) < 1e-9
    @test s.killing_integrated > 0
    @test s.detachment_integrated > 0
    @test s.clipped_biomass >= 0
    @test all(isfinite, s.phi_final)
    @test minimum(s.phi_final) >= 0
    @test length(s.z) == 65
end
