# SDDE Solvers

`solve(prob::AbstractSDDEProblem, alg; kwargs)`

Solves the SDDE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

The recommended method for SDDE problems are the `SDE` algorithms. On SDEs you
simply reuse the same algorithm as the `SDE` solver, and StochasticDelayDiffEq.jl
will convert it to an SDDE solver. The recommendations for SDDE solvers match
those of SDEs, except that only up to strong order 1 is recommended. This is
because the theoretical issues with higher order methods on SDDEs is currently
unknown. 
