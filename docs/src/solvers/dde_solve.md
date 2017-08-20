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

### Nonstiff DDEs

The standard choice is `MethodOfSteps(OrwenZen5())`. This is a highly efficient
FSAL 5th order algorithm which optimizes its interpolation error and should
handle most problems. For fast solving at where non-strict error control is
needed, choosing `OrwenZen3()` can do well. Using `BS3` is similar to the MATLAB
`dde23`, but `OrwenZen3()` will have noticably less error for the same work.
For algorithms where strict error control is needed, it is recommended that one
uses `DP8()`. Other high order integrators are not applicable since they use
a lazy interpolant.

For state-dependent delays, the current best choice is `RK4` since it uses a
residual control method to more accurately step. Note that current state-dependent
delays are not detected and thus non-residual control methods will be less
accurate. Still, residual control is an error-prone method. We recommend setting
the tolerances low (`1e-10`) and only trusting the solution to a 2-3 decimal
places of accuracy.

### Stiff DDEs and Differential-Algebraic Delay Equations (DADEs)

For stiff DDEs, the SDIRK and Rosenbrock methods are very efficient as they will
reuse the Jacobian in the unconstrained stepping iterations. One should choose
from the methods which have stiff-aware interpolants for better stability.
`MethodOfSteps(Rosenbrock23())` is a good low order method choice. Additionally,
the `Rodas` methods like `MethodOfSteps(Rodas4())` are good choices because of
their higher order stiff-aware interpolant.

Additionally, DADEs can be solved by specifying the problem in mass matrix form.
The Rosenbrock methods are good choices in these situations.

### Undeclared Lags

Lags are declared separately from their use. One can use any lag by simply using
the interpolant of `h` at that point. However, one should use caution in order
to achieve the best accuracy. When lags are declared, the solvers can more
efficiently be more accurate. If there are undeclared lags, one should only
use residual control methods like `RK4()` as these will better detect the
discontinuities.

## Special Keyword Arguments

- `minimal_solution` - Allows the algorithm to delete past history when `dense`
  and `save_everystep` are true. Defaults to true. If lags can grow this may
  need to be set to false.

### Note

If the method is having trouble, one may want to adjust the parameters of the
fixed-point iteration. Decreasing the absolute tolerance `fixedpoint_abstol` and the
relative tolerance `fixedpoint_reltol`, and increasing the maximal number of iterations
`max_fixedpoint_iters` can help ensure that the steps are correct. If the problem still
is not correctly converging, one should lower `dtmax`. In the worst case scenario, one
may need to set `constrained=true` which will constrain timesteps to at most the size
of the minimal lag and hence forces more stability at the cost of smaller timesteps.

There is currently no recommended algorithm for state-dependent delay problems.
An algorithm is currently in the works.
