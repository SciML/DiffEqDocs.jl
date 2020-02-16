# [Split ODE Problems](@id split_ode_prob)

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
\frac{du}{dt} =  Au + f_2(u,p,t)
```

For how to define a linear function `A`, see the documentation for the [DiffEqOperators](@ref).

### Constructors

```julia
SplitODEProblem(f::SplitFunction,u0,tspan,p=NullParameters();kwargs...)
SplitODEProblem{isinplace}(f1,f2,u0,tspan,p=NullParameters();kwargs...)
```

The `isinplace` parameter can be omitted and will be determined using the signature of `f2`.
Note that both `f1` and `f2` should support the in-place style if `isinplace` is `true` or they
should both support the out-of-place style if `isinplace` is `false`. You cannot mix up the two styles.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

Under the hood, a `SplitODEProblem` is just a regular `ODEProblem` whose `f` is a `SplitFunction`.
Therefore you can solve a `SplitODEProblem` using the same solvers for `ODEProblem`. For solvers
dedicated to split problems, see [Split ODE Solvers](@ref split_ode_solve).

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](@ref performance_overloads)
page.

### Fields

* `f1`, `f2`: The functions in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `p`: The parameters for the problem. Defaults to `NullParameters`
* `kwargs`: The keyword arguments passed onto the solves.
