# Chapter 7 — Optimal Bounds on Dynamical Systems via Auxiliary Functions

This directory contains the MATLAB scripts accompanying Chapter 7 of

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**  
> Jason Bramburger

The examples compute rigorous bounds on long-time averages, pointwise extrema, finite-time extreme events, maximal Lyapunov exponents, synchronization indicators, and controlled time-average quantities.

## Requirements

All scripts require:

- MATLAB
- YALMIP
- MOSEK

Some scripts also use standard MATLAB toolboxes such as Optimization Toolbox for nonlinear solves.

## Scripts

| Script | Description |
|---|---|
| `lorenz_time_average_bounds.m` | Computes SOS auxiliary-function bounds for normalized moments of the Lorenz system. |
| `henon_heiles_lyapunov_exponent.m` | Computes SOS bounds on the maximal Lyapunov exponent for the Hénon–Heiles system and generates the associated orbit figures. |
| `kuramoto_locking_bounds.m` | Computes auxiliary-function bounds for the Kuramoto locking residual. |
| `kuramoto_locked_branches.m` | Computes locked-state branches and bifurcation data for finite Kuramoto rings. |
| `KuramotoRing.m` | Shared helper class used by the Kuramoto scripts. |
| `ks_background.m` | Computes background-method and SOS bounds for the Kuramoto–Sivashinsky Galerkin model. |
| `rucklidge_pointwise_bounds.m` | Computes pointwise upper and lower bounds on the Rucklidge attractor. |
| `linear_extreme_event_bounds.m` | Computes finite-time extreme-event bounds for a nonnormal linear system. |
| `shear_flow_extreme_event_bounds.m` | Computes finite-time energy-growth bounds for the shifted MFE9 shear-flow model. |
| `controlled_vdp_duffing_bounds.m` | Synthesizes a controller to reduce an auxiliary-function bound for coupled Van der Pol–Duffing oscillators. |

## Notes

Some higher-degree SOS programs can be computationally expensive. The default parameters are chosen to make the scripts runnable while still reproducing the computations used in the book. Larger sweeps, higher polynomial degrees, and denser Monte Carlo simulations can be enabled by editing the options at the top of each script.

## Output

Scripts write figures, tables, and `.mat` files to the current working directory.
