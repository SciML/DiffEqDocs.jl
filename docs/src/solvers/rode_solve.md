# RODE Solvers

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package | Description |
|---|---|
| `StochasticDiffEqRODE` | Random ODE methods (`RandomEM`, `RandomTamedEM`, `RandomHeun`). |


## Recommended Methods

Currently, the only implemented method is the `RandomEM` method in StochasticDiffEq.jl.
It is strong order alpha for an alpha-Holder continuous noise process.

## Full List of Methods

### StochasticDiffEq.jl

Each of the StochasticDiffEq.jl solvers come with a linear interpolation.

  - `StochasticDiffEqRODE.RandomEM`- The Euler-Maruyama method for RODEs. Strong order matching Holder continuity.

Example usage:

```julia
using StochasticDiffEq    # RODEProblem, RandomEM
sol = solve(prob, RandomEM(), dt = 1 / 100)
```

!!! note "v8: load StochasticDiffEq directly"

    Under DifferentialEquations.jl v8 the umbrella only re-exports
    `OrdinaryDiffEq`, so `RandomEM` and the `RODEProblem` constructor must be
    obtained from `StochasticDiffEq` directly.
