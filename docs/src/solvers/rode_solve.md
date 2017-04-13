# RODE Solvers

## Recommended Methods

Currently, the only implemented method is the `RandomEM` method in StochasticDiffEq.jl.
It is strong order alpha for a alpha-Holder continuous noise process.

## Special Keyword Arguments

## Implemented Solvers

### StochasticDiffEq.jl

Each of the StochasticDiffEq.jl solvers come with a linear interpolation.

- `EM`- The Euler-Maruyama method for RODEs. Strong order matching Holder continuity.
