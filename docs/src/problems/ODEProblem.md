# Defining an ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE

```math
\frac{du}{dt} = f(t,u)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`),and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`, one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

Wraps the data which defines an ODE problem

```math
\\frac{du}{dt} = f(t,u)
```

with initial condition ``u0``.

### Constructors

`ODEProblem(f,u0,tspan)` : Defines the ODE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the ODE.
* `u0`: The initial condition.
* `isinplace`: Determines whether the function `f` uses the in-place syntax `f(t,u,du)`
  or not, `f(t,u)`
* `tspan`: The timespan for the problem.

## Example Problems

Examples problems can be found in [src/premades/premade_problems.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
prob = prob_ode_linear
sol = solve(prob,[0;1])
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
