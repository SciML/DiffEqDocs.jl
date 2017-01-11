# DDE Types

## Mathematical Specification of a DDE Problem

To define a DDE Problem, you simply need to give the function ``f`` and the initial
condition ``u0`` which define an ODE

```math
du = f(t,u,h)
```

`f` should be specified as `f(t,u,h)` (or in-place as `f(t,u,h,du)`).
`h` is the history function which is accessed for all delayed values. For example,
the `i`th component delayed by a time `tau` is denoted by `h(t-tau)`.
Note that we are not limited to numbers or vectors for `u0`, one is allowed to
provide `u0` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

```julia
ConstantLagDDEProblem(f,h,u0,lags,tspan)
DDEProblem(f,h,u0,lags,tspan)
```

### Fields

* `f`: The function in the ODE.
* `h`: The history function for the ODE before `t0`.
* `lags`: An array of lags. For constant lag problems this should be numbers.
  For state-dependent delay problems this is a tuple of functions.
* `tspan`: The timespan for the problem.

## Special Solver Options

## Special Solution Fields
