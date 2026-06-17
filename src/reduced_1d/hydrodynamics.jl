using Statistics

struct HydrodynamicConfig
    name::Symbol
    reynolds_min::Float64
    reynolds_max::Float64
    Da_values::Vector{Float64}
    Ktilde::Float64
    alpha_o::Float64
    L0_over_Lc::Float64
    Sc::Float64
end

function HydrodynamicConfig(name::Symbol)
    if name === :open
        return HydrodynamicConfig(:open, 1.0e3, 1.0e6, [60.0, 180.0, 600.0], 2.857142857142857, 1.111, 1.25e-3, 1000.0)
    elseif name === :micro
        return HydrodynamicConfig(:micro, 1.0e-1, 1.0e2, [6.0, 18.0, 60.0], 2.0, 3.333, 0.10, 1000.0)
    else
        error("Unknown hydrodynamic configuration: $(name). Use :open or :micro.")
    end
end

smoothstep(x::Real) = begin
    y = clamp(Float64(x), 0.0, 1.0)
    y*y*(3.0 - 2.0*y)
end

function log_blend_weight(Re::Real; Re_lo::Float64=3.0e4, Re_hi::Float64=3.0e5)
    return smoothstep((log10(Float64(Re)) - log10(Re_lo)) / (log10(Re_hi) - log10(Re_lo)))
end

sherwood_laminar(Re::Real, Sc::Real) = 0.664 * sqrt(Float64(Re)) * Float64(Sc)^(1.0/3.0)
sherwood_turbulent(Re::Real, Sc::Real) = 0.037 * Float64(Re)^(4.0/5.0) * Float64(Sc)^(1.0/3.0)

function sherwood_effective(config::HydrodynamicConfig, Re::Real)
    if config.name === :micro
        return sherwood_laminar(Re, config.Sc)
    end
    w = log_blend_weight(Re)
    return (1.0 - w) * sherwood_laminar(Re, config.Sc) + w * sherwood_turbulent(Re, config.Sc)
end

function biot_from_reynolds(config::HydrodynamicConfig, Re::Real)
    return config.alpha_o * config.L0_over_Lc * sherwood_effective(config, Re)
end

function skin_friction_laminar(Re::Real)
    return 1.328 / sqrt(Float64(Re))
end

function skin_friction_turbulent(Re::Real)
    return 0.0592 / Float64(Re)^(1.0/5.0)
end

function skin_friction_effective(Re::Real)
    w = log_blend_weight(Re)
    return (1.0 - w) * skin_friction_laminar(Re) + w * skin_friction_turbulent(Re)
end

function stress_proxy(config::HydrodynamicConfig, Re::Real, L::Real; χH::Float64=0.10, Φc_max::Float64=0.95)
    if config.name === :micro
        Φc = min(χH * Float64(L), Φc_max)
        return Float64(Re) * (1.0 - Φc)^(-2.0)
    end
    return 0.5 * Float64(Re)^2 * skin_friction_effective(Re)
end

function stress_reference(config::HydrodynamicConfig; samples::Int=300, Φc_ref::Float64=0.20)
    logs = range(log10(config.reynolds_min), log10(config.reynolds_max), length=samples)
    if config.name === :micro
        return maximum((10.0^x) * (1.0 - Φc_ref)^(-2.0) for x in logs)
    end
    return maximum(stress_proxy(config, 10.0^x, 1.0) for x in logs)
end

function detachment_coefficient(
    config::HydrodynamicConfig,
    Re::Real,
    L::Real;
    kd_max::Float64=0.013,
    q::Float64=0.5,
    σ_ref::Union{Nothing,Float64}=nothing,
    χH::Float64=0.10,
    Φc_max::Float64=0.95,
    Φc_ref::Float64=0.20,
)
    ref = σ_ref === nothing ? stress_reference(config; Φc_ref=Φc_ref) : σ_ref
    σ = stress_proxy(config, Re, L; χH=χH, Φc_max=Φc_max)
    return kd_max * (σ / ref)^q
end

function logspace(a::Real, b::Real, n::Integer)
    return collect(10.0 .^ range(log10(Float64(a)), log10(Float64(b)), length=Int(n)))
end
