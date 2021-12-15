# [ODE Problems](@id ode_prob)

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

- `ODEProblem(f::ODEFunction,u0,tspan,p=NullParameters();kwargs...)`
- `ODEProblem{isinplace}(f,u0,tspan,p=NullParameters();kwargs...)` :
  Defines the ODE with the specified functions. `isinplace` optionally sets whether
  the function is inplace or not. This is determined automatically, but not inferred.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](@ref performance_overloads)
page.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `p`: The parameters.
* `kwargs`: The keyword arguments passed onto the solves.

## Example Problems

Example problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/tree/master/src/ode).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.ODEProblemLibrary
# load problems
ODEProblemLibrary.importodeproblems()
prob = ODEProblemLibrary.prob_ode_linear
sol = solve(prob)
```

```@meta
CurrentModule = ODEProblemLibrary
```

```@docs
prob_ode_linear
prob_ode_2Dlinear
prob_ode_bigfloatlinear
prob_ode_bigfloat2Dlinear
prob_ode_large2Dlinear
prob_ode_2Dlinear_notinplace
prob_ode_lotkavoltera
prob_ode_fitzhughnagumo
prob_ode_threebody
prob_ode_pleiades
prob_ode_vanderpol
prob_ode_vanderpol_stiff
prob_ode_rober
prob_ode_rigidbody
prob_ode_hires
prob_ode_orego
prob_ode_pollution
prob_ode_nonlinchem
prob_ode_brusselator_1d
prob_ode_brusselator_2d
prob_ode_filament
prob_ode_thomas
prob_ode_lorenz
prob_ode_aizawa
prob_ode_dadras
prob_ode_chen
prob_ode_rossler
prob_ode_rabinovich_fabrikant
prob_ode_sprott
prob_ode_hindmarsh_rose
```
