# Discrete Problems

## Mathematical Specification of a Discrete Problem

To define an Discrete Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define a function map:

```math
u_{n+1} = f(u,p,t_{n+1})
```

`f` should be specified as `f(u,p,t)` (or in-place as `f(du,u,p,t)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well. `t_{n+1}` is the
current time at which the map is applied. For a `FunctionMap` with defaults,
`t_n = t0 + n*dt` (with `dt=1` being the default). For continuous-time Markov chains
this is the time at which the change is occuring.

Note that if the discrete solver is set to have `scale_by_time=true`, then the problem
is interpreted as the map:

```math
u_{n+1} = u_n + dt f(u,p,t_n)
```

## Problem Type

### Constructors

- `DiscreteProblem{isinplace}(f::ODEFunction,u0,tspan,p=NullParameters();kwargs...)` :
  Defines the discrete problem with the specified functions.
- `DiscreteProblem{isinplace}(f,u0,tspan,p=NullParameters();kwargs...)` :
  Defines the discrete problem with the specified functions.
- `DiscreteProblem{isinplace}(u0,tspan,p=NullParameters();kwargs...)` :
  Defines the discrete problem with the identity map.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](http://docs.juliadiffeq.org/latest/features/performance_overloads)
page.

### Fields

* `f`: The function in the map.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `p`: The parameters for the problem. Defaults to `NullParameters`
* `kwargs`: The keyword arguments passed onto the solves.

#### Note About Timing

Note that if no `dt` and not `tstops` is given, it's assumed that `dt=1` and thus
`tspan=(0,n)` will solve for `n` iterations. If in the solver `dt` is given, then
the number of iterations will change. And if `tstops` is not empty, the solver will
revert to the standard behavior of fixed timestep methods, which is "step to each
tstop".
