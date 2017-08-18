# DDE Solvers

`solve(prob::AbstractDDEProblem,alg;kwargs)`

Solves the DDE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

The recommended method for standard constant lag DDE problems are the `MethodOfSteps`
algorithms. These are constructed from an OrdinaryDiffEq.jl algorithm as follows:

```julia
MethodOfSteps(alg;constrained=false,
             fixedpoint_abstol = nothing,
             fixedpoint_reltol = nothing,
             fixedpoint_norm   = nothing,
             max_fixedpoint_iters = 10)
```

where `alg` is an OrdinaryDiffEq.jl algorithm. Most algorithms will work, though
a notable exception are algorithms which use a lazy interpolant (the Verner methods).

The standard choice is `MethodOfSteps(OrwenZen5())`. This is a highly efficient
FSAL 5th order algorithm which optimizes its interpolation error and should
handle most problems. For fast solving at where non-strict error control is
needed, choosing `OrwenZen3()` can do well. Using `BS3` is similar to the MATLAB
`dde23`, but `OrwenZen3()` will have noticably less error for the same work.
For algorithms where strict error control is needed, it is recommended that one
uses `DP8()`. Other high order integrators are not applicable since they use
a lazy interpolant.

If the method is having trouble, one may want to adjust the parameters of the
fixed-point iteration. Decreasing the absolute tolerance `fixedpoint_abstol` and the
relative tolerance `fixedpoint_reltol`, and increasing the maximal number of iterations
`max_fixedpoint_iters` can help ensure that the steps are correct. If the problem still
is not correctly converging, one should lower `dtmax`. In the worst case scenario, one
may need to set `constrained=true` which will constrain timesteps to at most the size
of the minimal lag and hence forces more stability at the cost of smaller timesteps.

There is currently no recommended algorithm for state-dependent delay problems.
An algorithm is currently in the works.
