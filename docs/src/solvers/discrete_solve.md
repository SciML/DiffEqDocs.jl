# Discrete Solvers

## DiscreteProblems

`solve(prob::DiscreteProblem,alg;kwargs)`

Solves the discrete function map defined by `prob` using the algorithm `alg`.
If no algorithm is given, a default algorithm will be chosen.

## Recommended Methods

The implementation for solving discrete equations is the `FunctionMap` algorithm
in OrdinaryDiffEq.jl. It allows the full common interface (including events/callbacks)
to solve function maps, along with everything else like plot recipes, while
completely ignoring the ODE functionality related to continuous equations (except
for a tiny bit of initialization). However, the `SimpleFunctionMap` from SimpleDiffEq.jl
can be more efficient if the mapping function is sufficiently cheap, but it doesn't have
all of the extras like callbacks and saving support (but does have an integrator interface).

## Full List of Methods

### OrdinaryDiffEq.jl

- `FunctionMap`: A basic function map which implements the full common interface.

OrdinaryDiffEq.jl also contains the `FunctionMap` algorithm which lets you 
It has a piecewise constant interpolation and allows for all of the 
callback/event handling capabilities (of course, with `rootfind=false`. If a 
`ContinuousCallback` is given, it's always assumed `rootfind=false`).

The constructor is:

```julia
FunctionMap()
FunctionMap{scale_by_time}()
```

Every step is the update

```math
u_{n+1} = f(t_{n+1},u_n).
```

If in addition `scale_by_time` is marked `true` (default is false), 
then every step is the update:

```math
u_{n+1} = u_n + dtf(t_{n+1},u_n).
```

Notice that this is the same as updates from the Euler method, except in this
case we assume that its a discrete change and thus the interpolation is
piecewise constant.

### SimpleDiffEq.jl

- `SimpleFunctionMap`: A barebones implementation of a function map. Is optimally-efficient
  and has an integrator interface version, but does not support callbacks or saving controls.
