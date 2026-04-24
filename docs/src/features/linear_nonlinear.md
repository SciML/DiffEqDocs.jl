# [Specifying (Non)Linear Solvers and Preconditioners](@id linear_nonlinear)

One of the key features of DifferentialEquations.jl is its flexibility. Keeping
with this trend, many of the native Julia solvers provided by DifferentialEquations.jl
allow you to choose the method for linear and nonlinear solving. This section
details how to make that choice.

!!! note
    
    We highly recommend looking at the [Solving Large Stiff Equations](@ref stiff)
    tutorial, which goes through these options in a real-world example.

!!! warning
    
    These options do not apply to the Sundials differential equation solvers
    (`CVODE_BDF`, `CVODE_Adams`, `ARKODE`, and `IDA`). For complete descriptions
    of similar functionality for Sundials, see the
    [Sundials ODE solver documentation](@ref ode_solve_sundials) and
    [Sundials DAE solver documentation](@ref dae_solve_sundials).

## Linear Solvers: `linsolve` Specification

For linear solvers, DifferentialEquations.jl uses
[LinearSolve.jl](https://github.com/SciML/LinearSolve.jl). Any
[LinearSolve.jl algorithm](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/)
can be used as the linear solver simply by passing the algorithm choice to
linsolve. For example, the following tells `TRBDF2` to use [KLU.jl](https://github.com/JuliaSparse/KLU.jl)

```julia
TRBDF2(linsolve = KLUFactorization())
```

Many choices exist, including GPU offloading, so consult the
[LinearSolve.jl documentation](https://docs.sciml.ai/LinearSolve/stable/) for more details
on the choices.

## Preconditioners: `precs` Specification

Any [LinearSolve.jl-compatible preconditioner](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/)
can be used as a left or right preconditioner. Starting with OrdinaryDiffEq
v7 / DiffEqBase v7, the per-solver `precs` kwarg and the old `DEFAULT_PRECS`
default were removed — preconditioners are now configured on the `linsolve`
object itself, through LinearSolve's `Pl` / `Pr` interface. For example:

```julia
using LinearSolve
alg = TRBDF2(linsolve = KrylovJL_GMRES(precs = mypreconditioner))
```

See the [LinearSolve preconditioners
documentation](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/)
for the expected `precs` function signature (`Pl, Pr = precs(A, p)`) and the
full list of supported preconditioner types.

## Nonlinear Solvers: `nlsolve` Specification

All the Julia-based implicit solvers (OrdinaryDiffEq.jl, StochasticDiffEq.jl, etc.)
allow for choosing the nonlinear solver that is used to handle the implicit system.
While fully modifiable and customizable, most users should stick to the pre-defined
nonlinear solver choices. These are:

  - `NLNewton(; κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5, new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true, relax = 0 // 1)`: A quasi-Newton method. The default.
  - `NLAnderson(; κ=1//100, max_iter=10, max_history::Int=5, aa_start::Int=1, droptol=nothing, fast_convergence_cutoff=1//5)`:
    Anderson acceleration. While more stable than functional iteration, this method
    is less stable than Newton's method, but does not require a Jacobian.
  - `NLFunctional(; κ=1//100, max_iter=10, fast_convergence_cutoff=1//5)`: This method
    is the least stable, but does not require a Jacobian. It should only be used for
    non-stiff ODEs.

The `NLNewton` solver allows for relaxation via the `relax` keyword parameter. Numerical values of `relax` must lie in the half open unit interval `[0,1)` where `0` corresponds to no relaxation. Alternatively, `relax` may be set to a line search algorithm from [LineSearches.jl](https://julianlsolvers.github.io/LineSearches.jl/stable/) in order to include a line search relaxation step in the Newton iterations. For example, the well known Newton-Armijo iterative scheme can be employed by setting `relax=BackTracking()` where `BackTracking` is provided by `LineSearches`.
