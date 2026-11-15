# RODE Solvers

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package                | Methods                                                     | Good for |
| ---------------------- | ----------------------------------------------------------- | -------- |
| `StochasticDiffEqRODE` | `RandomEM`, `RandomTamedEM`, `RandomHeun`, `RandomTaylor15` | Random ODEs (RODEs); time-dependent random forcing. |


## Recommended Methods

`RandomEM` is the default choice. On a RODE driven by Wiener noise it has strong
order 1 under the conditions below, and `RandomHeun` and `RandomTamedEM` measure the
same order there, so the three differ in cost and error constant rather than in
rate. `RandomHeun` evaluates `f` twice per step and is second order when `f` does
not depend on the noise. `RandomTamedEM` keeps every step bounded, for problems
where `RandomEM` blows up. When the noise is a stored path, a `NoiseGrid` resolved
more finely than the solver steps, `RandomTaylor15` reaches a higher order by using
the path between the step endpoints. See the [RODE tutorial](@ref rode_example) for
worked examples.

## Full List of Methods

### StochasticDiffEq.jl

Each of these solvers uses a fixed time step, so `dt` must be given, and comes with
a linear interpolation. If the problem sets no `noise`, a Wiener process starting at
0 drives it; any [noise process](@ref noise_process) can be passed instead, except
for `RandomTaylor15`, which needs a `NoiseGrid`.

  - `StochasticDiffEqRODE.RandomEM` - The Euler method for RODEs, with `W` taken at
    the start of each step. Strong order 1 under the conditions of Theorem 7.1 of
    Kloeden and Rosa, [arXiv:2306.15418](https://arxiv.org/abs/2306.15418), which
    cover Wiener and other semimartingale noise and require, among other things, `f`
    to be globally Lipschitz in `u` with a constant that does not depend on `t` or
    `W`, and to have bounded first and second derivatives in `W`. Fixed time step
    only.
  - `StochasticDiffEqRODE.RandomHeun` - A two-stage Heun method, with `W` taken at
    the start and at the end of each step. Measured strong order 1 on Wiener noise,
    and second order when `f` does not depend on `W`. Fixed time step only.
  - `StochasticDiffEqRODE.RandomTamedEM` - A tamed Euler method. Each step is
    ``u + dt\,k / (1 + dt\,\|k\|)`` with ``k = f(u, p, t, W)``, so each step moves
    `u` by at most about 1 in norm, however large `f` is. Measured strong order 1 on
    Wiener noise. Fixed time step only.
  - `StochasticDiffEqRODE.RandomTaylor15` - A derivative-free form of the order 1.5
    RODE-Taylor scheme of Han and Kloeden, *Random Ordinary Differential Equations
    and Their Numerical Solution* (Springer, 2017), Chapter 12, which has pathwise
    order 1.5 for Wiener noise, up to an arbitrarily small loss. The derivatives
    of `f` in `W` are replaced by finite differences with spacing `sqrt(dt)`, the
    drift is advanced by a Heun step, and the integrals of the path over each step
    are computed from the grid points inside it. Four evaluations of `f` per step.
    Measured strong order 1.7 to 2.0 on Wiener paths, and second order when `f`
    does not depend on `W`. Requires a `NoiseGrid` with scalar values. If no grid
    point lies inside a step, the order drops to 1 and a warning is given. Fixed
    time step only.

Example usage:

```julia
using StochasticDiffEq    # RODEProblem, RandomEM, RandomHeun, RandomTamedEM
sol = solve(prob, RandomEM(), dt = 1 / 100)
```

`RandomTaylor15` on a Wiener path stored on a grid finer than the solver steps:

```julia
using StochasticDiffEq    # RODEProblem, NoiseGrid, RandomTaylor15
t = collect(range(0.0, 1.0, length = 2^14 + 1))
W = [0.0; cumsum(sqrt(t[2]) .* randn(2^14))]
f(u, p, t, W) = -u * cos(W)
prob = RODEProblem(f, 1.0, (0.0, 1.0), noise = NoiseGrid(t, W))
sol = solve(prob, RandomTaylor15(), dt = 1 / 64)
```

!!! note "v8: load StochasticDiffEq directly"

    Under DifferentialEquations.jl v8 the umbrella only re-exports
    `OrdinaryDiffEq`, so `RandomEM`, `RandomHeun`, `RandomTamedEM`,
    `RandomTaylor15` and the `RODEProblem` constructor must be obtained from
    `StochasticDiffEq` directly.
