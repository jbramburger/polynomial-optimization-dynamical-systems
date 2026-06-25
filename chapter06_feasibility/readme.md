# Chapter 6 — Sum-of-Squares Certificates for Stability, Invariance, and Control

This directory contains the MATLAB scripts accompanying Chapter 6 of

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**  
> Jason Bramburger

The examples illustrate how sum-of-squares optimization can be used to certify stability, construct forward-invariant and trapping regions, approximate regions of attraction, exclude recurrent dynamics, and synthesize polynomial feedback controllers.

## Requirements

All scripts require

- MATLAB
- YALMIP
- MOSEK

Running each script reproduces the numerical results and figures appearing in the corresponding examples of the book.

## Scripts

| Script | Description |
|---|---|
| `sos_lyapunov_homogeneous_cubic_full.m` | Computes polynomial Lyapunov functions for a homogeneous cubic system using SOS optimization. |
| `sparse_lyapunov_scaling_demo.m` | Demonstrates the computational advantages of exploiting sparsity in SOS Lyapunov searches. |
| `backward_vdp_roa_sos.m` | Computes SOS inner approximations of the region of attraction for the backward Van der Pol oscillator. |
| `perturbed_hopf_3d_roa_sos.m` | Computes SOS inner approximations of the region of attraction for a perturbed three-dimensional Hopf system. |
| `lorenz_no_periodic_orbits.m` | Certifies parameter ranges for which the Lorenz equations admit no nontrivial periodic orbits inside a compact absorbing set. |
| `spiral_quadratic_trapping.m` | Computes certified quadratic trapping regions for a planar spiral system with an attracting periodic orbit. |
| `minimum_wave_speed.m` | Computes rigorous SOS upper and lower bounds on the minimum travelling-wave speed in a non-KPP isothermal diffusion model. |
| `pendulum_sos_controller_synthesis.m` | Synthesizes a polynomial feedback controller for the upright torque-actuated pendulum using an SOS Lyapunov certificate. |
| `vdp_clf_cbf_sos_synthesis.m` | Synthesizes polynomial feedback controllers for the controlled Van der Pol oscillator using combined CLF–CBF SOS certificates. |

## Notes

Some of the SOS programs become computationally intensive at higher polynomial degrees. The default parameters are chosen so that each script can be run on a typical desktop computer while reproducing the figures in the book. Polynomial degrees, solver tolerances, and problem parameters can be modified at the top of each script to explore stronger certificates or alternative examples.

## Output

Depending on the script, running an example produces

- printed numerical results,
- PDF figures,
- CSV tables of computed quantities, and
- MATLAB data files (where appropriate).

All output files are written to the current working directory.
