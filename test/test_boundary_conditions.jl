@testset "Boundary conditions" begin
    for mode in (:dirichlet, :robin)
        _, L, s = run_case(;
            N=33, dt=0.001, Tfinal=0.1,
            Lambda=0.0, Da=0.0, Da_a=0.0,
            C_init=0.7, A_init=0.12,
            C_top_bc=mode, C_top=0.7, C_top_robin_ref=0.7,
            C_top_robin_Bi=3.0,
            A_top_bc=mode, A_top=0.12, A_step_time=1.0,
            A_top_robin_Bi=3.0, write_output=false,
            fixed_thickness=true)
        @test L == 0.8
        @test maximum(abs.(s.C_final .- 0.7)) < 1e-11
        @test maximum(abs.(s.A_final .- 0.12)) < 1e-11
        @test maximum(abs.(s.phi_final .- 0.6)) < 1e-11
    end

    _, _, s = run_case(;
        N=41, Tfinal=0.1, dt=0.001, Lambda=0.0, Da=0.0,
        C_init=0.0, C_top_bc=:robin, C_top_robin_Bi=2.0,
        write_output=false, fixed_thickness=true)
    @test minimum(s.C_final) >= 0
    @test s.C_final[1] < s.C_final[end] < 1.0
end
