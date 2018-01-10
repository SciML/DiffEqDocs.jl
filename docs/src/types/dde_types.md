# DDE Problems

## Mathematical Specification of a DDE Problem

To define a DDE Problem, you simply need to give the function ``f`` and the initial
condition ``u0`` which define an ODE:

```math
du = f(t,u,h)
```

`f` should be specified as `f(t,u,h)` (or in-place as `f(t,u,h,du)`).
`h` is the history function which is accessed for all delayed values. For example,
the `i`th component delayed by a time `tau` is denoted by `h(t-tau)`.
Note that we are not limited to numbers or vectors for `u0`; one is allowed to
provide `u0` as arbitrary matrices / higher dimension tensors as well.

### Functional Forms for `h`

`h`, the history function, can be called the following ways:

- `h(t)`: out-of-place
- `h(out,t)` : in-place
- `h(t,deriv)` and `h(out,t,deriv)` where `deriv=Val{i}` is the in/out of place
  `i`th derivative calculation. To set this up, use `deriv::Type{Val{i}}` in the function
  signature.
- `h(t,deriv,idxs)` and `h(out,t,deriv,idxs)` where `idxs` is an integer for which
  index of the history to return.
  
Note that a dispatch for the supplied history function of matching form is required 
for whichever function forms are used in the user derivative function `f`. 

## Declaring Lags

Lags are declared separately from their use. One can use any lag by simply using
the interpolant of `h` at that point. However, one should use caution in order
to achieve the best accuracy. When lags are declared, the solvers can more
efficiently be more accurate and thus this is recommended.

## Problem Type

### Constructors

```julia
DDEProblem{isinplace}(f,h,u0,tspan,constant_lags=nothing,dependent_lags=nothing;
                      callback=nothing,mass_matrix=I)
```

`isinplace` optionally sets whether the function is inplace or not. This is
determined automatically, but not inferred.

### Fields

* `f`: The function in the ODE.
* `h`: The history function for the ODE before `t0`.
* `tspan`: The timespan for the problem.
* `constant_lags`: An array of constant lags. These should be numbers corresponding
  to times that are used in the history function `h`.
* `dependent_lags` A tuple of functions for the state-dependent lags used by the
  history function `h`.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.
