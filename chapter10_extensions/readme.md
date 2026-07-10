# Chapter 10: Extensions and Frontiers

This directory contains the MATLAB scripts accompanying the examples from Chapter 10 of

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**
>
> Jason Bramburger

## Requirements

All scripts in this directory require:

- MATLAB
- YALMIP
- MOSEK

Running each script reproduces the numerical examples and figures appearing in the corresponding chapter of the book.

## Scripts

| Script | Description |
| :--- | :--- |
| `duffing_escape_rate_bounds.m` | Computes upper and lower bounds on the principal eigenvalue governing escape from a stochastic Duffing double-well potential and reconstructs the leading eigenfunction from the SOS certificate. |
| `hybrid_clutch_control.m` | Synthesizes minimum-gain linear feedback controllers for a three-mode hybrid dry clutch model using SOS Lyapunov certificates, and compares controlled and uncontrolled switching dynamics. |
| `rational_sos_recovery.m` | Demonstrates rational recovery of an exact SOS Gram matrix from a numerical semidefinite solution and verifies the resulting certificate using exact arithmetic. |

## Notes

This chapter highlights several emerging directions for auxiliary-function methods beyond the classical settings developed earlier in the book. The examples include stochastic escape-rate estimation, Lyapunov-based controller synthesis for hybrid dynamical systems, and certified post-processing of numerical semidefinite programs through rational recovery. Together they illustrate how auxiliary-function techniques continue to expand into new application areas while bridging numerical optimization with mathematically rigorous certification.
