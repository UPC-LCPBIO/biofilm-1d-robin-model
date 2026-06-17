using DelimitedFiles
using Printf

# ============================================================
# Analytical asymptotic nutrient profiles for Figure 1
# Stationary first-order uptake limit: C << K
# ============================================================

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTDIR = joinpath(REPO_ROOT, "results", "raw", "result_1_nutrient_transfer")
mkpath(OUTDIR)

const outfile = joinpath(OUTDIR, "fig1_C_asymptotic_profiles.csv")

# -----------------------------
# Parameters from mechanism test
# -----------------------------
const Hb = 0.8
const Da = 180.0
const phi_l0 = 0.6
const Ktilde = 2.857142857142857

const Bi_values = [
    ("0p1", 0.1),
    ("1",   1.0),
    ("30",  30.0),
]

# -----------------------------
# First-order stationary profile
# -----------------------------
lambda_C = Da * phi_l0
m_C = sqrt(lambda_C / Ktilde)

function C_asymptotic(z_norm::Float64, Bi::Float64)
    Cs = Bi / (Bi + m_C * tanh(m_C * Hb))
    return Cs * cosh(m_C * Hb * z_norm) / cosh(m_C * Hb)
end

# -----------------------------
# Write CSV
# -----------------------------
z_norm = range(0.0, 1.0; length=501)

open(outfile, "w") do io
    println(io, "z_norm,C_asymp_Bi_0p1,C_asymp_Bi_1,C_asymp_Bi_30")

    for x in z_norm
        vals = [C_asymptotic(float(x), Bi) for (_, Bi) in Bi_values]
        @printf(io, "%.8f,%.10e,%.10e,%.10e\n", x, vals[1], vals[2], vals[3])
    end
end

println("Analytical profiles written to:")
println(outfile)
println()
println("Parameters:")
println("Hb      = ", Hb)
println("Da      = ", Da)
println("phi_l0  = ", phi_l0)
println("Ktilde  = ", Ktilde)
println("lambda_C = ", lambda_C)
println("m_C      = ", m_C)