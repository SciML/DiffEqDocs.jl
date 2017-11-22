# BVP Problems

## Mathematical Specification of an BVP Problem

To define an BVP Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
\frac{du}{dt} = f(t,u)
```

along with an implicit function `bc!` which defines the residual equation, where

```math
bc(t,u) = 0
```

is the manifold on which the solution must live. A common form for this is the
two-point `BVProblem` where the manifold defines the solution at two points:

```math
u(t_0) = a
u(t_f) = b
```

## Problem Type

### Constructors

```julia
TwoPointBVProblem{isinplace}(f,bc!,u0,tspan)
BVProblem{isinplace}(f,bc!,u0,tspan)
```

For `TwoPointBVProblem`, `bc!` is the inplace function:

```julia
bc!(residual, ua, ub)
```

where `residual` computed from the current ``u_a = u(t_0)`` and ``u_b = u(t_f)``.
For `BVProblem`, `bc!` is the inplace function:

```julia
bc!(residual, sol)
```

where `u` is the current solution to the ODE which is used to compute the `residual`.
Note that all features of the `ODESolution` are present in this form.
In both cases, the size of the residual matches the size of the initial condition.
For more general problems, use the
[parameter estimation routines](../../analysis/parameter_estimation.html) and change the signature of `bc!` and `f` to

```julia
bc!(residual, sol, p)
f(t, u, p)
```

where `p`is the `Vector` of free parameters.

### Fields

* `f`: The function for the ODE.
* `bc`: The boundary condition function.
* `u0`: The initial condition. Either the initial condition for the ODE as an
  initial value problem, or a `Vector` of values for ``u(t_i)`` for collocation
  methods
* `tspan`: The timespan for the problem.
* `p`: `Vector` of free parameters for parameter estimation
