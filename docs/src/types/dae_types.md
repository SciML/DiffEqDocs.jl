# DAE Problems

## Mathematical Specification of an DAE Problem

To define a DAE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
0 = f(du,u,p,t)
```

`f` should be specified as `f(du,u,p,t)` (or in-place as `f(resid,du,u,p,t)`).
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

- `DAEProblem(f::DAEFunction,du0,u0,tspan,p=NullParameters();kwargs...)`
- `DAEProblem{isinplace}(f,du0,u0,tspan,p=NullParameters();kwargs...)` :
  Defines the DAE with the specified functions.
  `isinplace` optionally sets whether the function is inplace or not. This is
  determined automatically, but not inferred.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](http://docs.juliadiffeq.org/latest/features/performance_overloads)
page.

### Fields

* `f`: The function in the ODE.
* `du0`: The initial condition for the derivative.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `differential_vars`: A logical array which declares which variables are the
  differential (non algebraic) vars (i.e. `du'` is in the equations for this
  variable). Defaults to nothing. Some solvers may require this be set if an
  initial condition needs to be determined.
* `p`: The parameters for the problem. Defaults to `NullParameters`
* `kwargs`: The keyword arguments passed onto the solves.

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/dae_premade_problems.jl).

To use a sample problem, such as `prob_dae_resrob`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.DAEProblemLibrary
# load problems
DAEProblemLibrary.importdaeproblems()
prob = DAEProblemLibrary.prob_dae_resrob
sol = solve(prob,IDA())
```

```@docs
DiffEqProblemLibrary.DAEProblemLibrary.prob_dae_resrob
```
