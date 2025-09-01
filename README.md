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
    \mathrm{minimize}   &   \quad   c^T x \\
    \mathrm{subject~to} &   \quad   Ax = b \\
                        &   \quad   x \geq 0,
\end{aligned}
$$

and [quadratic programs](https://en.wikipedia.org/wiki/Quadratic_programming), e.g.

$$
\begin{aligned}
    \mathrm{minimize}   &   \quad   \dfrac12 x^T P x + x^T q + r \\
    \mathrm{subject~to} &   \quad   Gx \leq h \\
                        &   \quad   Ax = b,
\end{aligned}
$$

where the inequalities are applied element-wise, and $P \in \mathbb{R}^{n \times n}$ is a symmetric positive semi-definite matrix.
While we strive for an intuitive user interface, we do not compromise with computational performances and make extensive use of Modern Fortran constructs for the implementation of battle-tested algorithms. The current version of `LightConvex` focuses on designing the structure of the package and its interface. For that purpose, `LightConvex` is thus currently restricted to linear programs.

#### Road map

### Building LightConvex

### Examples

### Licence

### Development

### References
