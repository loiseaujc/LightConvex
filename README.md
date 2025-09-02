# LightConvex

Lightweight package for [convex programming](https://en.wikipedia.org/wiki/Convex_optimization) written in Modern Fortran.

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/loiseaujc/LightConvex.svg)](https://github.com/loiseaujc/LightConvex/releases/latest)
[![Build Status](https://github.com/loiseaujc/LightConvex/actions/workflows/CI-CD.yml/badge.svg)](https://github.com/loiseaujc/LightConvex/actions)
[![codecov](https://codecov.io/gh/loiseaujc/LightConvex/branch/main/graph/badge.svg)](https://codecov.io/gh/loiseaujc/LightConvex)
[![last-commit](https://img.shields.io/github/last-commit/loiseaujc/LightConvex)](https://github.com/loiseaujc/LightConvex/commits/main)

**Status :** Prototype

## Description

`LightConvex` is an experiment about writing convex programming solvers in Modern Fortran. It aims at providing an easy-to-use API for solving both dense and sparse convex programs, including [linear programs](https://en.wikipedia.org/wiki/Linear_programming), e.g.

$$
\begin{aligned}
    \mathrm{maximize}   &   \quad   c^T x \\
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


## Road map

### Linear Programming

*Linear programming* (LP) is a method to achieve the best outcome (such as maximizing profit or minimizing cost) in a mathematical model whose requirements and objective are represented by linear relationships. The standard form of an LP reads

$$
\begin{aligned}
    \mathrm{maximize}   &   \quad   c^T x   \\
    \mathrm{subject~to} &   \quad   Ax = b  \\
                        &   \quad   x \geq 0,
\end{aligned}
$$

with $c \in \mathbb{R}^n$, $x \in \mathbb{R}^n$, $A \in \mathbb{R}^{m \times n}$ (with $m \leq n$) and $b \in \mathbb{R}^m$. Although it may seem restrictive at first, all LP can be formulated in this format.

Below is a tentative list of features to be included in the first working prototype of `LightConvex`:

- **Support for dense LP**
    - [x] Standard primal simplex algorithm
    - [ ] Standard dual simplex algorithm
    - [ ] Revised dual simplex algorithm
- **Support for sparse LP**
    - [ ] Primal simplex algorithm
    - [ ] Revised dual simplex algorithm
    - [ ] Primal Affine Scaling
- **High-level interfaces**
    - [ ] `problem = LP(c, A, b)` for dense and sparse LP.
    - [ ] `solution = solve(problem, alg)` for dense and sparse LP.
- **Preprocessing**
    - [ ] Conversion of any LP into standard equality form
    - [ ] Idiot Crash Algorithm
- **Utilities**
    - [ ] (Compressed) [MPS](https://en.wikipedia.org/wiki/MPS_(format)) file reader
- **Examples**
    - [ ] Max Flow / Min Cut
- **Documentation**
    - [ ] In-code documentation
    - [ ] Online documentation using `FORD`
- **Continuous integration and documentation**
    - [ ] Unit tests based on the [netlib LP test suite](https://www.netlib.org/lp)
    - [x] CI based on [setup-fortran-conda](https://github.com/gha3mi/setup-fortran-conda) with `fpm` build system and automatic documentation with `FORD`
    - [ ] Add code coverage using `codecov`

For the sake of simplicity, only double precision arithmetic is currently supported.

## Building LightConvex

### Fortran Package Manager (`fpm`)

The library can be build with the [Fortran Package Manager](https://github.com/fortran-lang/fpm) `fpm` using the provided `fpm.toml` like so:

```bash
    fpm build --release
```

Only double precision (`real64`) is currently supported.

To use `LightConvex` within your `fpm` project, add the following line to your `fpm.toml` file:

```toml
[dependencies]
LightConvex = {git="https://github.com/loiseaujc/LightConvex.git"}
```

## References

- Boyd, Stephen P., and Lieven Vandenberghe. Convex optimization. Cambridge university press, 2004.
