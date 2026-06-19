# Chapter 6: Sum-of-Squares Certificates for Stability, Invariance, and Control

This directory contains the MATLAB scripts accompanying the examples from Chapter 6 of

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**
>
> Jason Bramburger

## Requirements

All scripts in this directory require:

- MATLAB
- YALMIP
- MOSEK

Running each script should reproduce the numerical results and figures appearing in the corresponding examples of the book.

## Scripts

### `sos_lyapunov_homogeneous_cubic_full.m`

Computes SOS Lyapunov functions for a homogeneous cubic planar system.

**Output**

- Feasibility results for quadratic, quartic, and sextic Lyapunov searches.
- Several PDF figures showing Lyapunov functions, contours, trajectories, and decay.

---

### `sparse_lyapunov_scaling_demo.m`

Demonstrates the computational benefit of exploiting sparsity in SOS Lyapunov searches.

**Output**

- Structured Gram matrix sparsity pattern.
- Dense-versus-sparse scaling data.
- PDF figures showing sparsity and solve-time comparisons.

---

### `backward_vdp_roa_sos.m`

Computes SOS inner approximations of the region of attraction for the backward Van der Pol oscillator.

**Output**

- Certified area estimates for several polynomial degrees.
- Figure saved as `backward_vdp_roa_c1_variableV.pdf`.

---

### `perturbed_hopf_3d_roa_sos.m`

Computes SOS inner approximations of the region of attraction for a perturbed three-dimensional Hopf-type system.

**Output**

- Certified volume estimates for several polynomial degrees.
- Figure saved as `perturbed_hopf_3d_roa_c1.png`.

---

### `lorenz_no_periodic_orbits.m`

Certifies parameter ranges for which the Lorenz equations have no nontrivial periodic orbits inside a compact absorbing set.

**Output**

- Certified ranges of the Lorenz parameter rho.
- Figure saved as `lorenz_no_periodic_certification.pdf`.

---

### `spiral_quadratic_trapping.m`

Computes quadratic trapping regions for a planar spiral system with an attracting periodic orbit.

**Output**

- Certified quadratic trapping ellipses.
- Figure saved as `spiral_sos_quadratic_trapping_boundaries.pdf`.

---

### `minimum_wave_speed.m`

Computes SOS upper and lower bounds for the minimum traveling-wave speed in a non-KPP isothermal diffusion model.

**Output**

- Numerical lower and upper bounds as functions of the diffusion parameter.
- Figure saved as `nonkpp_minimum_wave_speed_bounds.pdf`.

---

### `pendulum_sos_controller_synthesis.m`

Synthesizes a polynomial feedback controller for the upright torque-actuated pendulum using an SOS Lyapunov certificate.

**Output**

- Synthesized polynomial feedback law.
- Figures saved as `pendulum_panelA.pdf` and `pendulum_panelB.pdf`.

---

### `vdp_clf_cbf_sos_synthesis.m`

Synthesizes polynomial feedback controllers for the controlled Van der Pol oscillator using combined CLF-CBF SOS certificates.

**Output**

- Synthesized controllers for two certified safe sets.
- Figures saved as `vdp_clf_cbf_ellipse.pdf` and `vdp_clf_cbf_lobe.pdf`.

## Notes

These examples illustrate SOS certificates for stability, trapping, region-of-attraction estimation, exclusion of recurrent dynamics, and feedback synthesis. The scripts are self-contained and are intended to be run directly from MATLAB. Readers are encouraged to modify polynomial degrees, solver tolerances, and system parameters to explore the behavior of the certificates.
