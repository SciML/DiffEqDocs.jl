# DDE Solvers

`solve(prob::AbstractDDEProblem,alg;kwargs)`

Solves the DDE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

The recommended method for standard constant lag DDE problems are the `MethodOfSteps`
algorithms. These are constructed from an OrdinaryDiffEq.jl algorithm as follows:

```julia
MethodOfSteps(alg;constrained=false,
             picardabstol = nothing,
             picardreltol = nothing,
             picardnorm   = nothing,
             max_picard_iters = 10)
```

where `alg` is an OrdinaryDiffEq.jl algorithm. Most algorithms will work, though
a notable exception are algorithms which use a lazy interpolant (the Verner methods).

The standard choice is `MethodOfSteps(Tsit5())`. This is a highly efficient FSAL
5th order algorithm which should handle most problems. For fast solving at where
non-strict error control is needed, choosing `BS3()` can do well (this is similar
to the MATLAB `dde23`). For algorithms where strict error control is needed, it
is recommended that one uses `DP8()`.

If the method is having trouble, one may want to adjust the Picard parameters.
Decreasing the Picard tolerances and increasing the Picard iterations can help
ensure that the steps are correct. If the problem still is not correctly converging,
one should lower `dtmax`. In the worst case scenarios, one may need to set
`constrained=true` which will constrain the method in a manner that forces
more stability at the cost of smaller timesteps.

There is currently no recommended algorithm for state-dependent delay problems.
An algorithm is currently in the works.
