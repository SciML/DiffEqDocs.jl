# SDDE Solvers

`solve(prob::AbstractSDDEProblem, alg; kwargs)`

Solves the SDDE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

The recommended method for SDDE problems are the `SDE` algorithms. On SDEs you
simply reuse the same algorithm as the `SDE` solver, and StochasticDelayDiffEq.jl
will convert it to an SDDE solver. The recommendations for SDDE solvers match
those of SDEs, except that only up to strong order 1 is recommended. Note too
that order 1 is currently only attainable if there is no delay term in the
diffusion function ``g``: delays in the drift function ``f`` are compatible
with first order convergence. Theoretical issues with higher order methods
(1.5+) on SDDEs is currently unknown.
