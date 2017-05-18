# Steady State Solvers

`solve(prob::SteadyStateProblem,alg;kwargs)`

Solves for the steady states in the problem defined by `prob` using the algorithm
`alg`. If no algorithm is given, a default algorithm will be chosen.

## Recommended Methods

Currently the only method is `SSRootfind` and so I am pretty sure it's the best
option right now.

# Full List of Methods

## SteadyStateDiffEq.jl

- `SSRootfind` : Using a rootfinding algorithm to find a steady state. Defaults
  to using NLsolve.jl. A different algorithm can be specified via the `nlsolve`
  keyword argument.
