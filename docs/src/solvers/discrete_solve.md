# Discrete Solvers

`solve(prob::DiscreteProblem,alg;kwargs)`

Solves the discrete function map defined by `prob` using the algorithm `alg`.
If no algorithm is given, a default algorithm will be chosen.

## Recommended Methods

The implementation for solving discrete equations is the `Discrete` algorithm
in OrdinaryDiffEq.jl. It has zero overhead and uses compilation to build a separate
setup that allows you to use the common interface (including events/callbacks)
to solve function maps, along with everything else like plot recipes, while
completely ignoring the ODE functionality related to continuous equations (except
for a tiny bit of initialization).

## Discrete Algorithm

OrdinaryDiffEq.jl also contains the `Discrete` algorithm which lets you solve
a problem where `f` is a map: ``u_{n+1} = f(t_{n+1},u_n)``. It has a piecewise constant
interpolation and allows for all of the callback/event handling capabilities
(of course, with `rootfind=false`. If a `ContinuousCallback` is given, it's always
assumed `rootfind=false`).

The constructor is:

```julia
Discrete(;apply_map=false,scale_by_time=false)
```

If `apply_map=false`, `f` is completely ignored. If `apply_map=true`, then
every step is the update

```math
u_{n+1} = f(t_{n+1},u_n).
```

If in addition `scale_by_time=true`,
then every step is the update

```math
u_{n+1} = u_n + dtf(t_{n+1},u_n).
```

Notice that this is the same as updates from the Euler method, except in this
case we assume that its a discrete change and thus the interpolation is piecewise constant.

As a shorthand,

```julia
FunctionMap(scale_by_time=false)
```

is a `Discrete` with `apply_map=true`, and thus corresponds to the function map
equation

```math
u_{n+1} = f(t_n,u_n).
```
