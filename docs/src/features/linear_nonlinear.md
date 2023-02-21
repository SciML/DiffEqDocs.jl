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
[LinearSolve.jl algorithm](https://linearsolve.sciml.ai/dev/solvers/solvers/)
can be used as the linear solver simply by passing the algorithm choice to
linsolve. For example, the following tells `TRBDF2` to use [KLU.jl](https://github.com/JuliaSparse/KLU.jl)

```julia
TRBDF2(linsolve = KLUFactorization())
```

Many choices exist, including GPU offloading, so consult the
[LinearSolve.jl documentation](https://linearsolve.sciml.ai/dev/) for more details
on the choices.

## Preconditioners: `precs` Specification

Any [LinearSolve.jl-compatible preconditioner](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/)
can be used as a left or right preconditioner. Preconditioners are specified by
the `Pl,Pr = precs(W,du,u,p,t,newW,Plprev,Prprev,solverdata)` function where
the arguments are defined as:

  - `W`: the current Jacobian of the nonlinear system. Specified as either
    ``I - \gamma J`` or ``I/\gamma - J`` depending on the algorithm. This will
    commonly be a `WOperator` type defined by OrdinaryDiffEq.jl. It is a lazy
    representation of the operator. Users can construct the W-matrix on demand
    by calling `convert(AbstractMatrix,W)` to receive an `AbstractMatrix` matching
    the `jac_prototype`.
  - `du`: the current ODE derivative
  - `u`: the current ODE state
  - `p`: the ODE parameters
  - `t`: the current ODE time
  - `newW`: a `Bool` which specifies whether the `W` matrix has been updated since
    the last call to `precs`. It is recommended that this is checked to only
    update the preconditioner when `newW == true`.
  - `Plprev`: the previous `Pl`.
  - `Prprev`: the previous `Pr`.
  - `solverdata`: Optional extra data the solvers can give to the `precs` function.
    Solver-dependent and subject to change.

The return is a tuple `(Pl,Pr)` of the LinearSolve.jl-compatible preconditioners.
To specify one-sided preconditioning, simply return `nothing` for the preconditioner
which is not used.

Additionally, `precs` must supply the dispatch:

```julia
Pl, Pr = precs(W, du, u, p, t, ::Nothing, ::Nothing, ::Nothing, solverdata)
```

which is used in the solver setup phase to construct the integrator
type with the preconditioners `(Pl,Pr)`.

The default is `precs=DEFAULT_PRECS` where the default preconditioner function
is defined as:

```julia
DEFAULT_PRECS(W, du, u, p, t, newW, Plprev, Prprev, solverdata) = nothing, nothing
```

## Nonlinear Solvers: `nlsolve` Specification

All the Julia-based implicit solvers (OrdinaryDiffEq.jl, StochasticDiffEq.jl, etc.)
allow for choosing the nonlinear solver that is used to handle the implicit system.
While fully modifiable and customizable, most users should stick to the pre-defined
nonlinear solver choices. These are:

  - `NLNewton(; κ=1//100, max_iter=10, fast_convergence_cutoff=1//5, new_W_dt_cutoff=1//5)`: A quasi-Newton method. The default.
  - `NLAnderson(; κ=1//100, max_iter=10, max_history::Int=5, aa_start::Int=1, droptol=nothing, fast_convergence_cutoff=1//5)`:
    Anderson acceleration. While more stable than functional iteration, this method
    is less stable than Newton's method, but does not require a Jacobian.
  - `NLFunctional(; κ=1//100, max_iter=10, fast_convergence_cutoff=1//5)`: This method
    is the least stable, but does not require a Jacobian. It should only be used for
    non-stiff ODEs.
