# Steady State Problems

## Mathematical Specification of a Steady State Problem

To define an Steady State Problem, you simply need to give the function ``f``
which defines the ODE:

```math
\frac{du}{dt} = f(u,p,t)
```

and an initial guess ``u₀`` of where `f(u,p,t)=0`. `f` should be specified as `f(u,p,t)`
(or in-place as `f(du,u,p,t)`), and `u₀` should be an AbstractArray (or number)
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
SteadyStateProblem(f::ODEFunction,u0)
SteadyStateProblem{isinplace}(f,u0)
```

`isinplace` optionally sets whether the function is inplace or not. This is
determined automatically, but not inferred. Additionally, the constructor from
`ODEProblem`s is provided:

```julia
SteadyStateProblem(prob::ODEProblem)
```

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](http://docs.juliadiffeq.org/latest/features/performance_overloads.html)
page.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial guess for the steady state.

## Special Solution Fields

The `SteadyStateSolution` type is different from the other DiffEq solutions because
it does not have temporal information.
