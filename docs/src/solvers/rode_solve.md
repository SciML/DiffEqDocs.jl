# RODE Solvers

## Recommended Methods

Currently, the only implemented method is the `RandomEM` method in StochasticDiffEq.jl.
It is strong order alpha for a alpha-Holder continuous noise process.

### StochasticDiffEq.jl

Each of the StochasticDiffEq.jl solvers come with a linear interpolation.

- `RandomEM`- The Euler-Maruyama method for RODEs. Strong order matching Holder continuity.
