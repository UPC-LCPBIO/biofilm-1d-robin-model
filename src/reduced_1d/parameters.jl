const L0tilde = 0.8
const N = 320
const Tfinal = 30.0
const dt = 2.0e-4
const save_every = 1000

const Ktilde = 2.857142857142857
const Lambda = 0.2
const Da = 180.0
const Da_a = 0.0
const beta_Ab = 0.2222222222222222
const PiE0 = 21.0
const kill_power = 1.0
const gamma_kill = 0.8
const live_frac0 = 0.6

const k_d = 0.0
const Lmin = 1.0e-4

const C_top_bc = :robin
const C_top = 1.0
const C_top_robin_Bi = 3.0
const C_top_robin_ref = 1.0

const A_top_bc = :robin
const A_top = 0.0
const A_top_robin_Bi = 3.0
const A_step_time = 6.0
const A_step_value = 0.041

const C_init = 0.0
const A_init = 0.0
const save_final_profiles = true

function threshold_case_kwargs(; Bi::Float64=3.0, A0::Float64=A_step_value)
    return (
        C_top_bc = :robin,
        C_top_robin_Bi = Bi,
        C_top_robin_ref = 1.0,
        A_top_bc = :robin,
        A_top_robin_Bi = Bi,
        A_step_time = 6.0,
        A_step_value = A0,
        PiE0 = 21.0,
    )
end

function mechanism_case_kwargs(; BiC::Float64=30.0, A0::Float64=0.11)
    return (
        C_top_bc = :robin,
        C_top_robin_Bi = BiC,
        C_top_robin_ref = 1.0,
        A_top_bc = :dirichlet,
        A_step_time = 6.0,
        A_step_value = A0,
        PiE0 = 54.0,
    )
end
