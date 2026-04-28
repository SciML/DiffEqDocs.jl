# DDE Solvers

`solve(prob::AbstractDDEProblem, alg; kwargs)`

Solves the DDE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package | Description |
|---|---|
| `DelayDiffEq` | DDE / SDDE driver: `MethodOfSteps`, `DDEProblem`, `SDDEProblem`. Wraps an inner ODE algorithm. |
| `OrdinaryDiffEqTsit5` | Tsit5 (5/4) and AutoTsit5 - the fastest non-stiff default for medium-tolerance problems. |
| `OrdinaryDiffEqVerner` | Vern6/7/8/9 and AutoVern - high-order explicit non-stiff methods with lazy variants. |
| `OrdinaryDiffEqLowOrderRK` | Low-order explicit Runge-Kutta methods (BS3, DP5, RK4, Heun, Euler, OwrenZen, etc.). |
| `OrdinaryDiffEqHighOrderRK` | High-order explicit RK beyond Verner (DP8, TanYam7, TsitPap8, PFRK87). |
| `OrdinaryDiffEqRosenbrock` | Rosenbrock and Rosenbrock-W methods for stiff problems (Rosenbrock23, Rodas4/5P, ROS variants). |
| `OrdinaryDiffEqSDIRK` | SDIRK / ESDIRK methods for stiff problems (KenCarp3/4/47/58, TRBDF2, ImplicitEuler, Kvaerno). |
| `OrdinaryDiffEqBDF` | BDF / NDF multistep (FBDF, QNDF, ABDF2, SBDF) and implicit-DAE forms (DFBDF, DImplicitEuler). |


## Recommended Methods

The recommended method for DDE problems are the `MethodOfSteps` algorithms.
These are constructed from an OrdinaryDiffEq.jl algorithm as follows:

```julia
MethodOfSteps(alg; constrained = false, fpsolve = NLFunctional(; max_iter = 10))
```

where `alg` is an OrdinaryDiffEq.jl algorithm. Most algorithms should work.

!!! note "v8: DelayDiffEq must be loaded explicitly"

    Under DifferentialEquations.jl v8 the `using DifferentialEquations`
    umbrella only re-exports `OrdinaryDiffEq`.  `MethodOfSteps`,
    `DDEProblem`, and the DDE-specific solver paths come from
    `DelayDiffEq.jl`; load them with

    ```julia
    using DelayDiffEq        # MethodOfSteps, DDEProblem
    # plus the OrdinaryDiffEq sublib(s) for the inner ODE algorithm, e.g.
    using OrdinaryDiffEqTsit5: Tsit5
    alg = MethodOfSteps(Tsit5())
    ```

    `DelayDiffEq` already `@reexport`s `OrdinaryDiffEq`'s default solver set
    (`Tsit5`, `Vern6`–`Vern9`, `Rosenbrock23`, `Rodas5P`, `FBDF`, ...), so
    `MethodOfSteps(Tsit5())` works after just `using DelayDiffEq`. Non-default
    inner algorithms (`BS3`, `RK4`, `DP8`, `Rodas4`, `Rodas5`, ...) require
    pulling in the corresponding `OrdinaryDiffEqXxx` sublib explicitly.

### Nonstiff DDEs

The standard algorithm choice is `MethodOfSteps(Tsit5())`. This is a highly efficient
FSAL 5th order algorithm with free interpolants which should
handle most problems. For fast solving where non-strict error control is
needed, choosing `MethodOfSteps(BS3())` can do well. Using `BS3` is similar to the MATLAB
`dde23`. For algorithms where strict error control is needed, it is recommended
that one uses `MethodOfSteps(Vern6())`. Benchmarks show that going to higher order methods like
`MethodOfSteps(DP8())` may not be beneficial.

### Stiff DDEs and Differential-Algebraic Delay Equations (DADEs)

For stiff DDEs, the SDIRK and Rosenbrock methods are very efficient, as they will
reuse the Jacobian in the unconstrained stepping iterations. One should choose
from the methods which have stiff-aware interpolants for better stability.
`MethodOfSteps(Rosenbrock23())` is a good low order method choice. Additionally,
the `Rodas` methods like `MethodOfSteps(Rodas4())` are good choices because of
their higher order stiff-aware interpolant.

Additionally, DADEs can be solved by specifying the problem in mass matrix form.
The Rosenbrock methods are good choices in these situations.

### Lag Handling

Lags are declared separately from their use. One can use any lag by simply using
the interpolant of `h` at that point. However, one should use caution in order
to achieve the best accuracy. When lags are declared, the solvers can be more
efficient and accurate. Constant delays are propagated until the
order is higher than the order of the integrator. If state-dependent delays are
declared, the algorithm will detect discontinuities arising from these delays and
adjust the step size such that these discontinuities are included in the mesh, if
steps are rejected. This way, all discontinuities are treated exactly.

If there are undeclared lags, the discontinuities due to delays are not tracked.
In this case, one should only use residual control methods like `MethodOfSteps(RK4())`,
which is the current best choice, as these will step more accurately.
Still, residual control is an error-prone method. We recommend setting the
tolerances lower in order to get accurate results, though this may be costly
since it will use a rejection-based approach to adapt to the delay discontinuities.

## Special Keyword Arguments

  - `discontinuity_interp_points` - Number of interpolation points used to track
    discontinuities arising from dependent delays. Defaults to 10. Only relevant
    if dependent delays are declared.

  - `discontinuity_abstol` and `discontinuity_reltol` - These are absolute and
    relative tolerances used by the check whether the time point at the beginning
    of the current step is a discontinuity arising from dependent delays. Defaults
    to 1/10^12 and 0. Only relevant if dependent delays are declared.

### Note

If the method is having trouble, one may want to adjust the fixed-point iteration.
Decreasing the absolute tolerance and the relative tolerance by specifying the
keyword arguments `abstol` and `reltol` when solving the DDE problem, and increasing
the maximal number of iterations by specifying the keyword argument `max_iter` in
the `MethodOfSteps` algorithm, can help ensure that the steps are correct. If the problem
still is not correctly converging, one should lower `dtmax`. For problems with only
constant delays, in the worst-case scenario, one may need to set `constrained = true` which
will constrain timesteps to at most the size of the minimal lag and hence forces more
stability at the cost of smaller timesteps.
