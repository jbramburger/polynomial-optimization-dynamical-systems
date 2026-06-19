# Polynomial Optimization Methods for Dynamical Systems

### Auxiliary Functions, Sum-of-Squares Programming, and Applications

This repository contains the MATLAB codes accompanying the book

> **Polynomial Optimization Methods for Dynamical Systems: Auxiliary Functions, Sum-of-Squares Programming, and Applications**
>
> Jason Bramburger

The repository includes scripts for reproducing the figures, computations, and numerical experiments presented throughout the text. The code is organized by chapter, with each directory containing the files associated with the corresponding examples and exercises.

## Software Requirements

Most examples require:

- MATLAB
- YALMIP
- MOSEK

Some examples additionally require:

- Chebfun

Specific requirements are indicated in the documentation for each chapter.

## Repository Structure

```text
chapter05_sos/           Polynomial optimization and SOS programming
chapter06_feasibility/   Lyapunov functions, regions of attraction,
                         trapping regions, and control certificates
chapter07_bounds/        Bounds on time averages, parameter-dependent
                         bounds, background methods, and extreme events
chapter08_measures/      Invariant measures, moment relaxations,
                         ergodic optimization, and duality
chapter09_data/          Data-driven polynomial optimization,
                         Koopman and Perron--Frobenius methods
chapter10_extensions/    Advanced topics and extensions
```

Each chapter directory contains a README file describing the associated
examples and listing the MATLAB scripts required to reproduce the figures
and computations appearing in the book.

## Citation

If you use the code in this repository, please cite:

```bibtex
@book{Bramburger2027,
  author    = {Jason Bramburger},
  title     = {Polynomial Optimization Methods for Dynamical Systems:
               Auxiliary Functions, Sum-of-Squares Programming, and Applications},
  publisher = {TBD},
  year      = {2027}
}
```

## License

The MATLAB code in this repository is released under the MIT License.

The accompanying book, its text, figures, and exercises remain subject to their respective copyright protections.
