# RODE Solvers

## Recommended Methods

Currently, the only implemented method is the `RandomEM` method in StochasticDiffEq.jl.
It is strong order alpha for an alpha-Holder continuous noise process.

## Full List of Methods

### StochasticDiffEq.jl

Each of the StochasticDiffEq.jl solvers come with a linear interpolation.

- `RandomEM`- The Euler-Maruyama method for RODEs. Strong order matching Holder continuity.

Example usage:

```julia
sol = solve(prob,RandomEM(),dt=1/100)
```
