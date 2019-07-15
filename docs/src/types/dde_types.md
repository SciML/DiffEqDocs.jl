# DDE Problems

## Mathematical Specification of a DDE Problem

To define a DDE Problem, you simply need to give the function ``f``, the initial
condition ``u_0`` at time point ``t_0``, and the history function ``h``
which together define a DDE:

```math
\begin{align*}
    \frac{du}{dt} &= f(u,h,p,t) \qquad & (t \geq t_0), \\
    u(t_0) &= u_0, \\
    u(t) &= h(t) \qquad &(t < t_0).
\end{align*}
```

``f`` should be specified as `f(u, h, p, t)` (or in-place as `f(du, u, h, p, t)`),
``u_0`` should be an AbstractArray (or number) whose geometry matches the
desired geometry of `u`, and ``h`` should be specified as described below. The
history function `h` is accessed for all delayed values. Note that we are not
limited to numbers or vectors for ``u_0``; one is allowed to provide ``u_0``
as arbitrary matrices / higher dimension tensors as well.

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
DDEProblem(f[, u0], h, tspan[, p]; <keyword arguments>)
DDEProblem{isinplace}(f[, u0], h, tspan[, p]; <keyword arguments>)
```

Parameter `isinplace` optionally sets whether the function is inplace or not.
This is determined automatically, but not inferred.

For specifying Jacobians and mass matrices, see the [DiffEqFunctions](http://docs.juliadiffeq.org/latest/features/performance_overloads.html) page.

### Arguments

* `f`: The function in the DDE.
* `u0`: The initial condition. Defaults to the value `h(p, first(tspan))` of the history function evaluated at the initial time point.
* `h`: The history function for the DDE before `t0`.
* `tspan`: The timespan for the problem.
* `p`: The parameters with which function `f` is called. Defaults to `nothing`.
* `constant_lags`: A collection of constant lags used by the history function `h`. Defaults to `()`.
* `dependent_lags` A tuple of functions `(u, p, t) -> lag` for the state-dependent lags
  used by the history function `h`. Defaults to `()`.
* `neutral`: If the DDE is neutral, i.e., if delays appear in derivative terms.
* `order_discontinuity_t0`: The order of the discontinuity at the initial time point. Defaults to `0` if an initial condition `u0` is provided. Otherwise it is forced to be greater or equal than `1`.
* `callback`: A callback to be applied to every solver which uses the problem. Defaults to `nothing`.
