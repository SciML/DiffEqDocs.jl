# RODE Solvers

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package                | Methods                                   | Good for |
| ---------------------- | ----------------------------------------- | -------- |
| `StochasticDiffEqRODE` | `RandomEM`, `RandomTamedEM`, `RandomHeun` | Random ODEs (RODEs); time-dependent random forcing. |


## Recommended Methods

`RandomEM` is the default choice. On a RODE driven by Wiener noise it has strong
order 1 under the conditions below, and `RandomHeun` and `RandomTamedEM` measure the
same order there, so the three differ in cost and error constant rather than in
rate. `RandomHeun` evaluates `f` twice per step and is second order when `f` does
not depend on the noise. `RandomTamedEM` keeps every step bounded, for problems
where `RandomEM` blows up. See the [RODE tutorial](@ref rode_example) for worked
examples.

## Full List of Methods

### StochasticDiffEq.jl

Each of these solvers uses a fixed time step, so `dt` must be given, and comes with
a linear interpolation. If the problem sets no `noise`, a Wiener process starting at
0 drives it; any [noise process](@ref noise_process) can be passed instead.

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

Example usage:

```julia
using StochasticDiffEq    # RODEProblem, RandomEM, RandomHeun, RandomTamedEM
sol = solve(prob, RandomEM(), dt = 1 / 100)
```

!!! note "v8: load StochasticDiffEq directly"

    Under DifferentialEquations.jl v8 the umbrella only re-exports
    `OrdinaryDiffEq`, so `RandomEM`, `RandomHeun`, `RandomTamedEM` and the
    `RODEProblem` constructor must be obtained from `StochasticDiffEq` directly.
