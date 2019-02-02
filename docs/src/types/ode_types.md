# ODE Problems

## Mathematical Specification of an ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
\frac{du}{dt} = f(u,p,t)
```

`f` should be specified as `f(u,p,t)` (or in-place as `f(du,u,p,t)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

- `ODEProblem(f::ODEFunction,u0,tspan,callback=CallbackSet())`
- `ODEProblem{isinplace}(f,u0,tspan,callback=CallbackSet())` :
  Defines the ODE with the specified functions. `isinplace` optionally sets whether
  the function is inplace or not. This is determined automatically, but not inferred.

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](http://docs.juliadiffeq.org/latest/features/performance_overloads.html)
page.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.

## Example Problems

Example problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/tree/master/src/ode).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary
prob = prob_ode_linear
sol = solve(prob)
```

```@docs
DiffEqProblemLibrary.prob_ode_linear
DiffEqProblemLibrary.prob_ode_2Dlinear
DiffEqProblemLibrary.prob_ode_bigfloatlinear
DiffEqProblemLibrary.prob_ode_bigfloat2Dlinear
DiffEqProblemLibrary.prob_ode_large2Dlinear
DiffEqProblemLibrary.prob_ode_2Dlinear_notinplace
DiffEqProblemLibrary.prob_ode_threebody
DiffEqProblemLibrary.prob_ode_pleides
DiffEqProblemLibrary.prob_ode_vanderpol
DiffEqProblemLibrary.prob_ode_vanderpol_stiff
DiffEqProblemLibrary.prob_ode_rober
DiffEqProblemLibrary.prob_ode_rigidbody
```
