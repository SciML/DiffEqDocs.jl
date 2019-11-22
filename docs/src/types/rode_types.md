# RODE Problems

## Mathematical Specification of a RODE Problem

To define a RODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
\frac{du}{dt} = f(u,p,t,W(t))
```

where `W(t)` is a random process. `f` should be specified as `f(u,p,t,W)`
(or in-place as `f(du,u,p,t,W)`), and `u₀` should be an AbstractArray (or number)
whose geometry matches the desired geometry of `u`. Note that we are not limited
to numbers or vectors for `u₀`; one is allowed to provide `u₀` as arbitrary matrices
/ higher dimension tensors as well.

### Constructors

- `RODEProblem(f::RODEFunction,u0,tspan,p=nothing;noise=WHITE_NOISE,rand_prototype=nothing,callback=nothing)`
- `RODEProblem{isinplace}(f,u0,tspan,p=nothing;noise=WHITE_NOISE,rand_prototype=nothing,callback=nothing,mass_matrix=I)` :
  Defines the RODE with the specified functions. The default noise is `WHITE_NOISE`.
  `isinplace` optionally sets whether the function is inplace or not. This is
  determined automatically, but not inferred.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](http://docs.juliadiffeq.org/dev/features/performance_overloads)
page.

### Fields

* `f`: The drift function in the SDE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `p`: The optional parameters for the problem. Defaults to `NullParameters`.
* `noise`: The noise process applied to the noise upon generation. Defaults to
  Gaussian white noise. For information on defining different noise processes,
  see [the noise process documentation page](../../../features/noise_process)
* `rand_prototype`: A prototype type instance for the noise vector. It defaults
  to `nothing`, which means the problem should be interpreted as having a noise
  vector whose size matches `u0`.
* `kwargs`: The keyword arguments passed onto the solves.
