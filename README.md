# Biofilm 1D Robin model

This repository contains the one-dimensional biofilm model used to generate the computational results of the manuscript on transport-controlled biofilm survival under flow.

The model describes nutrient transport, antibiotic transport, live-biomass dynamics, growth-induced velocity and, when enabled, stress-based detachment through a moving biofilm thickness. The interfacial exchange with the surrounding fluid is represented through Robin-type transfer conditions.

## Repository structure

```text
biofilm-1d-robin-model/
├── Project.toml
├── README.md
├── scripts/
│   ├── reproduce_paper.jl
│   ├── result_1_nutrient_transfer.jl
│   ├── result_1_asymptotic_profiles.jl
│   ├── result_2_transfer_regime_map.jl
│   ├── result_3_reynolds_threshold.jl
│   └── result_4_stress_detachment.jl
├── src/
│   ├── BiofilmRobinModel.jl
│   ├── parameters.jl
│   ├── grid.jl
│   ├── linear_solvers.jl
│   ├── transport.jl
│   ├── growth_velocity.jl
│   ├── hydrodynamics.jl
│   ├── solver.jl
│   └── campaigns.jl
└── test/
    ├── runtests.jl
    ├── test_conservation.jl
    ├── test_boundary_conditions.jl
    ├── test_moving_boundary.jl
    ├── test_regression.jl
    └── test_convergence.jl
```

The package is loaded as `using BiofilmRobinModel`. Numerical integration is
implemented once in `src/solver.jl`; the paper scripts configure the model and
write experiment-specific diagnostics.

The `results/` directory is intentionally empty in the repository. The scripts create the required output folders when they are executed.

## Model summary

The biofilm is represented on a one-dimensional domain

$$
0 \leq \tilde z \leq \tilde L(\tilde t)
$$

where $\tilde C$ is the nutrient concentration, $\tilde A$ is the antibiotic concentration, $\phi_\ell$ is the live-biomass fraction, $\tilde v$ is the growth-induced velocity and $\tilde L$ is the biofilm thickness.

The nondimensional governing equations are

```math
\frac{\partial \tilde C}{\partial \tilde t}
=
\frac{\partial^2 \tilde C}{\partial \tilde z^2}
-
Da\,e(\tilde C)\,\phi_\ell
```

```math
\frac{\partial \tilde A}{\partial \tilde t}
=
\beta_a\frac{\partial^2 \tilde A}{\partial \tilde z^2}
-
D_{aa}\tilde A\phi_\ell
```

```math
\frac{\partial \phi_\ell}{\partial \tilde t}
+
\frac{\partial}{\partial \tilde z}\left(\tilde v\phi_\ell\right)
=
\Lambda e(\tilde C)\phi_\ell
-
\Lambda\Pi_E(e)\tilde A\phi_\ell^{\gamma_{\mathrm{kill}}}
```

with

```math
\frac{\partial \tilde v}{\partial \tilde z}
=
\Lambda e(\tilde C)\phi_\ell,
\qquad
\frac{d\tilde L}{d\tilde t}
=
\tilde v(\tilde L,\tilde t)-k_d\tilde L^2
```

The nutrient activation and antibiotic killing modulation are

$$
e(\tilde C)=\frac{\tilde C}{\tilde C+\tilde K},
\qquad
\Pi_E(e)=\Pi_{E0}e^p
$$

At the substratum, zero-flux conditions are imposed for nutrient and antibiotic, and the biomass velocity is set to zero. At the biofilm-fluid interface, nutrient and antibiotic exchange are imposed through Robin-type transfer laws when the corresponding transfer mode is enabled.

## Numerical method

The moving physical domain is mapped to the fixed computational interval using
$\xi=\tilde z/\tilde L(\tilde t)$. The live-biomass balance is written in
conservative arbitrary Lagrangian–Eulerian (ALE) form:

$$
\frac{\partial (\tilde L\phi_\ell)}{\partial\tilde t}
+\frac{\partial}{\partial\xi}
\left[(\tilde v-\xi\dot{\tilde L})\phi_\ell\right]
=\tilde L\left(\Lambda e(\tilde C)\phi_\ell
-\Lambda\Pi_{E0}e(\tilde C)^p\tilde A\phi_\ell^{\gamma_{\mathrm{kill}}}\right).
$$

Node-centered finite volumes use half-width boundary control volumes. Shared
upwind fluxes give a telescoping discrete biomass balance. At the moving
surface the outward relative velocity is the detachment speed
$\tilde v(\tilde L)-\dot{\tilde L}=k_d\tilde L^2$.

Nutrient and antibiotic equations use the corresponding ALE transformation,
with implicit diffusion, explicit geometric fluxes and reactions, and Robin
exchange implemented as a boundary flux in the top half-volume. Dirichlet
boundaries are supported as well. The time integrator is first-order IMEX;
upwind biomass transport is first-order in space. No remapping between physical
grids is needed.

For `fixed_thickness=true`, the original nonexpanding, nonadvecting biomass
model is retained for the first three paper studies. Because that reduced
configuration suppresses volumetric expansion, it does **not** generally enforce
$\phi_\ell\leq1$ under unrestricted growth; this is a modeling limitation,
not a guarantee of the numerical method. The moving-boundary discretization
should also be monitored for small violations of local volume-fraction bounds.

Biomass diagnostics returned by `run_case` include integrated growth, killing,
detachment and any amount restored by positivity clipping. The conservation
residual accounts for that last correction; a small residual alone does not
establish positivity, convergence or physical validity.

### Automated checks

Run the package tests with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests cover discrete conservation, moving-mesh geometric consistency,
boundary conditions, the output interface and diffusion-grid refinement.
For research figures, perform additional timestep and grid-refinement checks
for the specific regimes being reported.

## Computational results

The repository contains one Julia script per main result.

### Result 1: nutrient-transfer mechanism

```bash
julia --project=. scripts/result_1_nutrient_transfer.jl
julia --project=. scripts/result_1_asymptotic_profiles.jl
```

This fixed-thickness study varies the nutrient Biot number while keeping the antibiotic forcing fixed. It generates time-series and profile CSV files under

```text
results/raw/result_1_nutrient_transfer/
```

### Result 2: $(Bi,\tilde A_0)$ regime map

```bash
julia --project=. scripts/result_2_transfer_regime_map.jl
```

This fixed-thickness sweep computes the post-treatment response index

$$
\mathcal R=\frac{\Phi(\tilde T)}{\Phi(\tilde t_{\mathrm{start}})}
$$

over the external-transfer and antibiotic-dose parameter space. The output is written to

```text
results/raw/result_2_regime_map/
```

### Result 3: Reynolds-dependent antibiotic thresholds

```bash
julia --project=. scripts/result_3_reynolds_threshold.jl
```

This fixed-thickness study maps Reynolds-dependent transfer into critical antibiotic thresholds. It performs repeated bracket-and-bisection searches for the neutral condition $\mathcal R=1$, so the full sweep can be computationally expensive.

The output is written to

```text
results/raw/result_3_reynolds_threshold/
```

### Result 4: stress-based detachment response

```bash
julia --project=. scripts/result_4_stress_detachment.jl
```

This study activates the moving-boundary equation and the stress-based detachment closure. The output is written to

```text
results/raw/result_4_stress_detachment/
```

## Full reproduction command

To run all Julia result scripts in sequence:

```bash
julia --project=. scripts/reproduce_paper.jl
```

The full workflow can take a long time, mainly because Result 3 performs many threshold searches. For routine checks, it is preferable to run individual result scripts.

## Fixed-thickness versus moving-boundary configuration

The intended configuration is:

```text
Result 1 -> fixed thickness
Result 2 -> fixed thickness
Result 3 -> fixed thickness
Result 4 -> moving boundary with stress-based detachment
```

This separation keeps the transport-survival mechanisms isolated in the first three studies and introduces growth-detachment coupling only in the final Reynolds-dependent detachment analysis.

## Requirements

The code is written in Julia and uses only standard-library dependencies listed in `Project.toml`.

A Julia version compatible with the project can be used through:

```bash
julia --project=.
```