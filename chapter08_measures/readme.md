# Chapter 8 — Invariant Measures, Ergodic Optimization, and Duality

This directory contains the MATLAB scripts accompanying Chapter 8 of

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**  
> Jason Bramburger

The examples illustrate how moment relaxations and occupation measures can be used to approximate invariant measures, compute optimal long-time averages, reconstruct physical measures, and transport probability measures under nonlinear dynamics.

## Requirements

All scripts require

- MATLAB
- YALMIP
- MOSEK

Running each script reproduces the numerical results and figures appearing in the corresponding examples of the book.

## Scripts

| Script | Description |
|---|---|
| `vdp_moment_sos_time_average_bounds.m` | Compares moment-SDP bounds and SOS auxiliary-function bounds for long-time averages of the Van der Pol oscillator. |
| `rossler_physical_measure_moment_sdp.m` | Approximates the physical invariant measure of the Rössler attractor using a moment SDP and reconstructs a polynomial density. |
| `lorenz_occupation_measure_transport.m` | Computes finite-time terminal measures for the Lorenz system using occupation-measure moment constraints and reconstructs terminal densities. |

## Notes

These examples illustrate optimization over invariant, physical, occupation, and terminal measures using semidefinite programming. The scripts are self-contained and are intended to be run directly from MATLAB. Higher moment degrees and longer simulation horizons may improve accuracy but can substantially increase runtime.

## Output

Depending on the script, running an example produces

- printed numerical results,
- PDF or PNG figures,
- CSV tables of computed quantities, and
- MATLAB data files (where appropriate).

All output files are written to the current working directory.
