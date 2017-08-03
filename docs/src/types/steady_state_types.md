# Steady State Problems

## Mathematical Specification of a Steady State Problem

To define an Steady State Problem, you simply need to give the function ``f``
which defines the ODE:

```math
\frac{du}{dt} = f(t,u)
```

and an initial guess ``u₀`` of where `f(t,u)=0`. `f` should be specified as `f(t,u)`
(or in-place as `f(t,u,du)`), and `u₀` should be an AbstractArray (or number)
whose geometry matches the desired geometry of `u`. Note that we are not limited
to numbers or vectors for `u₀`; one is allowed to provide `u₀` as arbitrary
matrices / higher dimension tensors as well.

Note that for the steady-state to be defined, we must have that `f` is autonomous,
that is `f` is independent of `t`. But the form which matches the standard ODE
solver should still be used. The steady state solvers interpret the `f` by
fixing `t=0`.

## Problem Type

### Constructors

```julia
SteadyStateProblem{isinplace}(f,u0,mass_matrix=I)
```

`isinplace` optionally sets whether the function is inplace or not. This is
determined automatically, but not inferred. Additionally, the constructor from
the `ODEProblem` is provided:

```julia
SteadyStateProblem(prob::ODEProblem)
```

### Fields

* `f`: The function in the ODE.
* `u0`: The initial guess for the steady state.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

## Special Solution Fields

The `SteadyStateSolution` type is different from the other DiffEq solutions because
it does not have temporal information.
