# Split ODE Problems

## Mathematical Specification of a Split ODE Problem

To define a `SplitODEProblem`, you simply need to give a two functions
functions ``f_1`` and ``f_2`` along with an initial condition ``u₀`` which
define an ODE:

```math
\frac{du}{dt} =  f_1(u,p,t) + f_2(u,p,t)
```

`f` should be specified as `f(u,p,t)` (or in-place as `f(du,u,p,t)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

Many splits are at least partially linear. That is the equation:

```math
\frac{du}{dt} =  f_1(u,p,t) + f_2(u,p,t)
```

For how to define a function as linear, see the documentation for the [DiffEqOperators](../../features/diffeq_operator.html).

### Constructors

```julia
SplitODEProblem{isinplace}(f1,f2,u0,tspan;kwargs...)
```

### Fields

* `f1`, `f2`: The functions in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.
