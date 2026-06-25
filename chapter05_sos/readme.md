# Chapter 5 — Polynomial Optimization and Sum-of-Squares Programming

This directory contains the MATLAB scripts accompanying Chapter 5 of

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**  
> Jason Bramburger

The examples introduce fundamental concepts in polynomial optimization and sum-of-squares programming, including global polynomial optimization, Putinar certificates on semialgebraic sets, and certification of sublevel set containment.

## Requirements

All scripts require

- MATLAB
- YALMIP
- MOSEK

Running each script reproduces the numerical results and figures appearing in the corresponding examples of the book.

## Scripts

| Script | Description |
|---|---|
| `motzkin_global_vs_putinar.m` | Illustrates the distinction between global SOS certificates and Putinar certificates using the Motzkin polynomial. Computes the best global SOS lower bound, the best Putinar lower bound on the unit disk, and visualizes the corresponding minimizers. |
| `sublevel_containment_example.m` | Computes the largest certified polynomial sublevel set contained within another semialgebraic set using a Putinar certificate and visualizes the certified containment. |

## Notes

These examples are intended to illustrate the construction and interpretation of sum-of-squares certificates rather than large-scale optimization. The scripts are self-contained and are designed to be run directly from MATLAB. Readers are encouraged to modify the polynomial functions, semialgebraic sets, and relaxation degrees to explore the behavior of the resulting certificates.

## Output

Depending on the script, running an example produces

- printed numerical results,
- PDF figures, and
- MATLAB data files (where appropriate).

All output files are written to the current working directory.
