# DAE Types

## Mathematical Specification of an DAE Problem

To define a DAE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
0 = f(t,u,du)
```

`f` should be specified as `f(t,u,du)` (or in-place as `f(t,u,du,resid)`).
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

`DAEProblem(f,u0,du0,tspan)` : Defines the ODE with the specified functions.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `du0`: The initial condition for the derivative.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to a black CallbackSet, which will have no effect.
  
## Special Solver Options

## Special Solution Fields

* `du`: The saved derivative values.

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/dae_premade_problems.jl).

To use a sample problem, such as `prob_dae_resrob`, you can do something like:

```julia
#Pkg.add("DiffEqProblemLibrary")
using DiffEqProblemLibrary
prob = prob_dae_resrob
sol = solve(prob,IDA())
```
