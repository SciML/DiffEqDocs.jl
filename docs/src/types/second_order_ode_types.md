# Second Order ODE Types

## Mathematical Specification of an Second Order ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
u'' = f(t,u,u')
```

`f` should be specified as `f(t,u,du)` (or in-place as `f(t,u,du,ddu)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

`SecondOrderODEProblem(f,u0,du0,tspan,callback=CallbackSet(),mass_matrix=I)` : Defines the ODE with the specified functions.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `du0`: The initial derivative.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

## Special Solver Options

## Special Solution Fields
