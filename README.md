# Biofilm 1D Robin reduced model

This repository contains the reduced one-dimensional biofilm model used to generate the computational results of the manuscript on transport-controlled biofilm survival under flow.

The model describes nutrient transport, antibiotic transport, live-biomass dynamics, growth-induced velocity and, when enabled, stress-based detachment through a moving biofilm thickness. The interfacial exchange with the surrounding fluid is represented through Robin-type transfer conditions.

## Repository structure

```text
biofilm-1d-robin-model/
├── README.md
├── LICENSE
├── Project.toml
├── .gitignore
├── results/
├── scripts/
│   ├── reproduce_paper.jl
│   ├── result_1_nutrient_transfer.jl
│   ├── result_1_asymptotic_profiles.jl
│   ├── result_2_transfer_regime_map.jl
│   ├── result_3_reynolds_threshold.jl
│   └── result_4_stress_detachment.jl
└── src/
    ├── BiofilmRobinReducedModel.jl
    └── reduced_1d/
        ├── ReducedBiofilm1D.jl
        ├── campaigns.jl
        ├── grid.jl
        ├── growth_velocity.jl
        ├── hydrodynamics.jl
        ├── linear_solvers.jl
        ├── main_1d_paper.jl
        ├── parameters.jl
        ├── solver.jl
        └── transport.jl
```

The `results/` directory is intentionally empty in the repository. The scripts create the required output folders when they are executed.

## Model summary

The biofilm is represented on a one-dimensional domain

\[
0 \leq \tilde z \leq \tilde L(\tilde t),
\]

where \(\tilde C\) is the nutrient concentration, \(\tilde A\) is the antibiotic concentration, \(\phi_\ell\) is the live-biomass fraction, \(\tilde v\) is the growth-induced velocity and \(\tilde L\) is the biofilm thickness.

The nondimensional governing equations are

\[
\frac{\partial \tilde C}{\partial \tilde t}
=
\frac{\partial^2 \tilde C}{\partial \tilde z^2}
-
Da\,e(\tilde C)\,\phi_\ell,
\]

\[
\frac{\partial \tilde A}{\partial \tilde t}
=
\beta_a\frac{\partial^2 \tilde A}{\partial \tilde z^2}
-
D_{aa}\tilde A\phi_\ell,
\]

\[
\frac{\partial \phi_\ell}{\partial \tilde t}
+
\frac{\partial}{\partial \tilde z}\left(\tilde v\phi_\ell\right)
=
\Lambda e(\tilde C)\phi_\ell
-
\Lambda\Pi_E(e)\tilde A\phi_\ell^{\gamma_{\mathrm{kill}}},
\]

with

\[
\frac{\partial \tilde v}{\partial \tilde z}
=
\Lambda e(\tilde C)\phi_\ell,
\qquad
\frac{d\tilde L}{d\tilde t}
=
\tilde v(\tilde L,\tilde t)-k_d\tilde L^2.
\]

The nutrient activation and antibiotic killing modulation are

\[
e(\tilde C)=\frac{\tilde C}{\tilde C+\tilde K},
\qquad
\Pi_E(e)=\Pi_{E0}e^p.
\]

At the substratum, zero-flux conditions are imposed for nutrient and antibiotic and the biomass velocity satisfies \(\tilde v(0,\tilde t)=0\). At the biofilm-fluid interface, nutrient and antibiotic exchange are imposed through Robin-type transfer laws when the corresponding transfer mode is enabled.

## Numerical method

The one-dimensional equations are solved on a uniform grid over the current biofilm domain. The implementation uses a first-order IMEX finite-difference strategy: diffusive transport is advanced implicitly through tridiagonal linear systems, while reaction terms and biomass advection are treated explicitly. The biomass advection term is discretised using an upwind flux.

The fixed-thickness studies impose

\[
\tilde L(\tilde t)=\tilde L(0),
\]

whereas the stress-detachment study activates the moving-boundary equation.

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

### Result 2: \((Bi,\tilde A_0)\) regime map

```bash
julia --project=. scripts/result_2_transfer_regime_map.jl
```

This fixed-thickness sweep computes the post-treatment response index

\[
\mathcal R=\frac{\Phi(\tilde T)}{\Phi(\tilde t_{\mathrm{start}})}
\]

over the external-transfer and antibiotic-dose parameter space. The output is written to

```text
results/raw/result_2_regime_map/
```

### Result 3: Reynolds-dependent antibiotic thresholds

```bash
julia --project=. scripts/result_3_reynolds_threshold.jl
```

This fixed-thickness study maps Reynolds-dependent transfer into critical antibiotic thresholds. It performs repeated bracket-and-bisection searches for the neutral condition \(\mathcal R=1\), so the full sweep can be computationally expensive.

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

