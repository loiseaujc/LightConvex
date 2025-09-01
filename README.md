# LightConvex

---

Lightweight package for [convex programming](https://en.wikipedia.org/wiki/Convex_optimization) written in Modern Fortran.

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
<!-- [![GitHub release](https://img.shields.io/github/release/loiseaujc/LightConvex.svg)](https://github.com/loiseaujc/LightConvex/releases/latest) -->
<!-- [![Build Status](https://github.com/loiseaujc/LightConvex/actions/workflows/ci.yml/badge.svg)](https://github.com/loiseaujc/LightConvex/actions) -->
<!-- [![codecov](https://codecov.io/gh/loiseaujc/LightConvex/branch/main/graph/badge.svg)](https://codecov.io/gh/loiseaujc/LightConvex) -->
[![last-commit](https://img.shields.io/github/last-commit/loiseaujc/LightConvex)](https://github.com/loiseaujc/LightConvex/commits/main)

**Status :** Prototype

### Description

`LightConvex` is an experiment about writing convex programming solvers in Modern Fortran. It aims at providing an easy-to-use API for solving both dense and sparse convex programs, including [linear programs](https://en.wikipedia.org/wiki/Linear_programming), e.g.

$$
\begin{aligned}
    \mathrm{minimize}   &   \quad   \mathbf{c}^T \mathbf{x} \\
    \mathrm{subject~to} &   \quad   \mathbf{Ax} = \mathbf{b} \\
                        &   \quad   \mathbf{x} \geq \mathbf{0},
\end{aligned}
$$

and [quadratic programs](https://en.wikipedia.org/wiki/Quadratic_programming), e.g.

$$
\begin{aligned}
    \mathrm{minimize}   &   \quad   \dfrac12 \mathbf{x}^\mathbf{P}x + \mathbf{x}^T \mathbf{q} + r \\
    \mathrm{subject~to} &   \quad   \mathbf{Gx} \leq \mathbf{h} \\
                        &   \quad   \mathbf{Ax} = \mathbf{b}.
\end{aligned}
$$

While we strive for an intuitive user interface, we do not compromise with computational performances and make extensive use of Modern Fortran constructs for the implementation of battle-tested algorithms. The current version of `LightConvex` focuses on designing the structure of the package and its interface. For that purpose, `LightConvex` is thus currently restricted to linear programs.

#### Road map

### Building LightConvex

### Examples

### Licence

### Development

### References
