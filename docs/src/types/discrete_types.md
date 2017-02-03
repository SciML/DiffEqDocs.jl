# Discrete Types

## Mathematical Specification of a Discrete Problem

To define an Discrete Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define a function map:

```math
u_{n+1} = f(t,u_n)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

Note that if the discrete solver is set to have `scale_by_time=true`, then the problem
is interpreted as the map:

```math
u_{n+1} = u_n + dtf(t,u_n)
```

## Problem Type

### Constructors

`DiscreteProblem(f,u0,tspan)` : Defines the discrete problem with the specified functions.

### Fields

* `f`: The function in the map.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.

## Special Solver Options

## Special Solution Fields

None. The Discrete type is as basic as it gets.
