# [SDDE Problems](@id sdde_prob)

## Mathematical Specification of a Stochastic Delay Differential Equation (SDDE) Problem

To define a SDDE Problem, you simply need to give the drift function ``f``,
the diffusion function `g`, the initial condition ``u_0`` at time point ``t_0``,
and the history function ``h`` which together define a SDDE:

```math
du = f(u,h,p,t)dt + g(u,h,p,t)dW_t \qquad & (t \geq t_0)
```
```math
u(t_0) = u_0,
```
```math
u(t) = h(t) \qquad &(t < t_0).
```

``f`` should be specified as `f(u, h, p, t)` (or in-place as `f(du, u, h, p, t)`)
(and ``g`` should match). ``u_0`` should be an AbstractArray (or number) whose
geometry matches the desired geometry of `u`, and ``h`` should be specified as
described below. The history function `h` is accessed for all delayed values.
Note that we are not limited to numbers or vectors for ``u_0``; one is allowed
to provide ``u_0`` as arbitrary matrices / higher dimension tensors as well.

Note that this functionality should be considered experimental.

## Functional Forms of the History Function

The history function `h` can be called in the following ways:

- `h(p, t)`: out-of-place calculation
- `h(out, p, t)`: in-place calculation
- `h(p, t, deriv::Type{Val{i}})`: out-of-place calculation of the `i`th derivative
- `h(out, p, t, deriv::Type{Val{i}})`: in-place calculation of the `i`th derivative
- `h(args...; idxs)`: calculation of `h(args...)` for indices `idxs`

Note that a dispatch for the supplied history function of matching form is required
for whichever function forms are used in the user derivative function `f`.

## Declaring Lags

Lags are declared separately from their use. One can use any lag by simply using
the interpolant of `h` at that point. However, one should use caution in order
to achieve the best accuracy. When lags are declared, the solvers can more
efficiently be more accurate and thus this is recommended.

## Problem Type

### Constructors

```
SDDEProblem(f,g[, u0], h, tspan[, p]; <keyword arguments>)
SDDEProblem{isinplace}(f,g[, u0], h, tspan[, p]; <keyword arguments>)
```

Parameter `isinplace` optionally sets whether the function is inplace or not.
This is determined automatically, but not inferred.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the [DiffEqFunctions](@ref performance_overloads) page.

### Arguments

* `f`: The drift function in the SDDE.
* `g`: The diffusion function in the SDDE.
* `u0`: The initial condition. Defaults to the value `h(p, first(tspan))` of the history function evaluated at the initial time point.
* `h`: The history function for the DDE before `t0`.
* `tspan`: The timespan for the problem.
* `p`: The parameters with which function `f` is called. Defaults to `NullParameters`.
* `constant_lags`: A collection of constant lags used by the history function `h`. Defaults to `()`.
* `dependent_lags` A tuple of functions `(u, p, t) -> lag` for the state-dependent lags
  used by the history function `h`. Defaults to `()`.
* `neutral`: If the DDE is neutral, i.e., if delays appear in derivative terms.
* `order_discontinuity_t0`: The order of the discontinuity at the initial time
  point. Defaults to `0` if an initial condition `u0` is provided. Otherwise
  it is forced to be greater or equal than `1`.
* `kwargs`: The keyword arguments passed onto the solves.
