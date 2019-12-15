# [DDE Problems](@id dde_prob)

## Mathematical Specification of a DDE Problem

To define a DDE Problem, you simply need to give the function ``f``, the initial
condition ``u_0`` at time point ``t_0``, and the history function ``h``
which together define a DDE:

```math
\frac{du}{dt} = f(u,h,p,t) \qquad & (t \geq t_0)
```
```math
u(t_0) = u_0,
```
```math
u(t) = h(t) \qquad &(t < t_0).
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

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the [DiffEqFunctions](@ref performance_overloads) page.

### Arguments

* `f`: The function in the DDE.
* `u0`: The initial condition. Defaults to the value `h(p, first(tspan))` of the history function evaluated at the initial time point.
* `h`: The history function for the DDE before `t0`.
* `tspan`: The timespan for the problem.
* `p`: The parameters with which function `f` is called. Defaults to `NullParameters`.
* `constant_lags`: A collection of constant lags used by the history function `h`. Defaults to `()`.
* `dependent_lags` A tuple of functions `(u, p, t) -> lag` for the state-dependent lags
  used by the history function `h`. Defaults to `()`.
* `neutral`: If the DDE is neutral, i.e., if delays appear in derivative terms.
* `order_discontinuity_t0`: The order of the discontinuity at the initial time point. Defaults to `0` if an initial condition `u0` is provided. Otherwise it is forced to be greater or equal than `1`.
* `kwargs`: The keyword arguments passed onto the solves.

## Example Problems

Example problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/tree/master/src/dde).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.ODEProblemLibrary
# load problems
ODEProblemLibrary.importodeproblems()
prob = ODEProblemLibrary.prob_ode_linear
sol = solve(prob)
```

### DDEs with 1 constant delay

```@meta
CurrentModule = DDEProblemLibrary
```

```@docs
prob_dde_constant_1delay_ip
prob_dde_constant_1delay_oop
prob_dde_constant_1delay_scalar
prob_dde_constant_1delay_long_ip
prob_dde_constant_1delay_long_oop
prob_dde_constant_1delay_long_scalar
```

### DDEs with 2 constant delays

```@docs
prob_dde_constant_2delays_ip
prob_dde_constant_2delays_oop
prob_dde_constant_2delays_scalar
prob_dde_constant_2delays_long_ip
prob_dde_constant_2delays_long_oop
prob_dde_constant_2delays_long_scalar
```

### DDETest Problems

Some details:

```
# DDEs with time dependent delays
prob_dde_DDETST_A1, prob_dde_DDETST_A2,
# DDEs with vanishing time dependent delays
prob_dde_DDETST_B1, prob_dde_DDETST_B2,
# DDEs with state dependent delays
prob_dde_DDETST_C1, prob_dde_DDETST_C2, prob_dde_DDETST_C3, prob_dde_DDETST_C4,
# DDEs with vanishing state dependent delays
prob_dde_DDETST_D1, prob_dde_DDETST_D2,
# neutral DDEs with time dependent delays
prob_dde_DDETST_E1, prob_dde_DDETST_E2,
# neutral DDEs with vanishing time dependent delays
prob_dde_DDETST_F1, prob_dde_DDETST_F2, prob_dde_DDETST_F3, prob_dde_DDETST_F4, prob_dde_DDETST_F5,
# neutral DDEs with state dependent delays
prob_dde_DDETST_G1, prob_dde_DDETST_G2,
# neutral DDEs with vanishing state dependent delays
prob_dde_DDETST_H1, prob_dde_DDETST_H2, prob_dde_DDETST_H3, prob_dde_DDETST_H4
```

```@docs
prob_dde_DDETST_A1
prob_dde_DDETST_A2
prob_dde_DDETST_B1
prob_dde_DDETST_B2
prob_dde_DDETST_C1
prob_dde_DDETST_C2
prob_dde_DDETST_C3
prob_dde_DDETST_C4
prob_dde_DDETST_D1
prob_dde_DDETST_D2
prob_dde_DDETST_E1
prob_dde_DDETST_E2
prob_dde_DDETST_F1
prob_dde_DDETST_F2
prob_dde_DDETST_F3
prob_dde_DDETST_F4
prob_dde_DDETST_F5
prob_dde_DDETST_G1
prob_dde_DDETST_G2
prob_dde_DDETST_H1
prob_dde_DDETST_H2
prob_dde_DDETST_H3
prob_dde_DDETST_H4
```

### Radar5 Test Problems

```@docs
prob_dde_RADAR5_oregonator
prob_dde_RADAR5_robertson
prob_dde_RADAR5_waltman
```

### QS Example

```@docs
prob_dde_qs
```
