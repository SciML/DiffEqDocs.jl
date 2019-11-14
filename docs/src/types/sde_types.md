# SDE Problems

## Mathematical Specification of a SDE Problem

To define an SDE Problem, you simply need to give the forcing function `f`,
the noise function `g`, and the initial condition `u₀` which define an SDE:

```math
du = f(u,p,t)dt + Σgᵢ(u,p,t)dWⁱ
```

`f` and `g` should be specified as `f(u,p,t)` and  `g(u,p,t)` respectively, and `u₀`
should be an AbstractArray whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well. A vector
of `g`s can also be defined to determine an SDE of higher Ito dimension.

## Problem Type

Wraps the data which defines an SDE problem

```math
u = f(u,p,t)dt + Σgᵢ(u,p,t)dWⁱ
```

with initial condition `u0`.

### Constructors

- `SDEProblem(f::SDEFunction,g,u0,tspan,p=nothing;noise=WHITE_NOISE,noise_rate_prototype=nothing)`
- `SDEProblem{isinplace}(f,g,u0,tspan,p=nothing;noise=WHITE_NOISE,noise_rate_prototype=nothing)` :
  Defines the SDE with the specified functions. The default noise is `WHITE_NOISE`.
  `isinplace` optionally sets whether the function is inplace or not. This is
  determined automatically, but not inferred.

Parameters are optional, and if not given then a `NullParameters()` singleton
will be used which will throw nice errors if you try to index non-existent
parameters. Any extra keyword arguments are passed on to the solvers. For example,
if you set a `callback` in the problem, then that `callback` will be added in
every solve call.

For specifying Jacobians and mass matrices, see the
[DiffEqFunctions](http://docs.juliadiffeq.org/latest/features/performance_overloads)
page.

### Fields

* `f`: The drift function in the SDE.
* `g`: The noise function in the SDE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `p`: The optional parameters for the problem. Defaults to `NullParameters`.
* `noise`: The noise process applied to the noise upon generation. Defaults to
  Gaussian white noise. For information on defining different noise processes,
  see [the noise process documentation page](../../../features/noise_process)
* `noise_rate_prototype`: A prototype type instance for the noise rates, that
  is the output `g`. It can be any type which overloads `A_mul_B!` with itself
  being the middle argument. Commonly, this is a matrix or sparse matrix. If
  this is not given, it defaults to `nothing`, which means the problem should
  be interpreted as having diagonal noise.  
* `kwargs`: The keyword arguments passed onto the solves.

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/sde_premade_problems.jl).

To use a sample problem, such as `prob_sde_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.SDEProblemLibrary
# load problems
SDEProblemLibrary.importsdeproblems()
prob = SDEProblemLibrary.prob_sde_linear
sol = solve(prob)
```

```@docs
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_linear
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_2Dlinear
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_wave
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_lorenz
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_cubic
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_additive
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_additivesystem
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_nltest
DiffEqProblemLibrary.SDEProblemLibrary.oval2ModelExample
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_stiffquadstrat
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_stiffquadito
DiffEqProblemLibrary.SDEProblemLibrary.generate_stiff_stoch_heat
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_bistable
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_bruss
DiffEqProblemLibrary.SDEProblemLibrary.prob_sde_oscilreact
```
