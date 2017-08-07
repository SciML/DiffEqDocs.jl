# Split ODE Problems

## Mathematical Specification of a Split ODE Problem

To define a `SplitODEProblem`, you simply need to give a tuple of
functions ``(f_1,f_2,\ldots,f_n)`` and the initial condition ``u₀`` which
define an ODE:

```math
\frac{du}{dt} =  f_1(t,u) + f_2(t,u) + \ldots + f_n(t,u)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

Many splits are at least partially linear. For how to define a function as linear,
see the documentation for the [DiffEqOperators](../../features/diffeq_operator.html).

### Constructors

```julia
SplitODEProblem(f1,...,fn,u0,tspan;kwargs...)
```

### Fields

* `f`: The functions in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.
