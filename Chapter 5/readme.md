# Chapter 5: Polynomial Optimization and Sum-of-Squares Programming

This directory contains the MATLAB scripts accompanying the examples from Chapter 5 of

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

### `motzkin_global_vs_putinar.m`

Illustrates the distinction between global SOS certificates and Putinar certificates on semialgebraic sets using the Motzkin polynomial

\[
m(x,y)=x^4y^2+x^2y^4+1-3x^2y^2.
\]

The script:

- Tests whether the Motzkin polynomial is SOS.
- Computes the best global SOS lower bound.
- Computes a Putinar lower bound on the unit disk.
- Produces a contour plot showing the global and constrained minimizers.

**Output**

- Numerical lower bounds printed to the MATLAB command window.
- Figure saved as `motzkin.pdf`.

---

### `sublevel_containment_example.m`

Demonstrates SOS certification of sublevel set containment.

Given two polynomial functions \(p(x,y)\) and \(q(x,y)\), the script computes the largest certified value \(c\) such that

\[
\{(x,y):p(x,y)\le c\}
\subseteq
\{(x,y):q(x,y)\le 1\}
\]

within a prescribed ambient semialgebraic set.

The script:

- Constructs a Putinar certificate for sublevel set containment.
- Uses bisection to determine the largest certified value of \(c\).
- Produces a visualization of the certified containment.

**Output**

- Certified value of \(c\).
- Figure saved as `sublevel_containment.pdf`.

## Notes

These examples are intended to illustrate the construction and interpretation of sum-of-squares certificates rather than large-scale optimization. Readers are encouraged to modify the polynomial functions, semialgebraic sets, and relaxation degrees to explore the behavior of the resulting certificates.
