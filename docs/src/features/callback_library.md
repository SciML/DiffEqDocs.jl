# Callback Library

DiffEqCallbacks.jl provides a library of various helpful callbacks which
can be used with any component solver which implements the callback interface.
It adds the following callbacks which are available to users of DifferentialEquations.jl.

## Manifold Conservation and Projection

In many cases, you may want to declare a manifold on which a solution lives.
Mathematically, a manifold `M` is defined by a function `g` as the set of points
where `g(u)=0`. An embedded manifold can be a lower dimensional object which
constrains the solution. For example, `g(u)=E(u)-C` where `E` is the energy
of the system in state `u`, meaning that the energy must be constant (energy
preservation). Thus by defining the manifold the solution should live on, you
can retain desired properties of the solution.

It is a consequence of convergence proofs both in the deterministic and stochastic
cases that post-step projection to manifolds keep the same convergence rate
(stochastic requires a truncation in the proof, details details), thus any algorithm
can be easily extended to conserve properties. If the solution is supposed to live
on a specific manifold or conserve such property, this guarantees the conservation
law without modifying the convergence properties.

### Constructor

```julia
ManifoldProjection(g; nlsolve=NLSOLVEJL_SETUP(), save=true, autonomous=numargs(g)==2, nlopts=Dict{Symbol,Any}())
```

- `g`: The residual function for the manifold. This is an inplace function of form
  `g(u, resid)` or `g(t, u, resid)` which writes to the residual the difference from
  the manifold components.
- `nlsolve`: A nonlinear solver as defined [in the nlsolve format](@ref linear_nonlinear).
- `save`: Whether to do the save after the callback is applied. Standard saving is unchanged.
- `autonomous`: Whether `g` is an autonomous function of the form `g(u, resid)`.
- `nlopts`: Optional arguments to nonlinear solver which can be any of the [NLsolve keywords](https://github.com/JuliaNLSolvers/NLsolve.jl#fine-tunings).

### Example

Here we solve the harmonic oscillator:

```julia
u0 = ones(2)
function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -u[1]
end
prob = ODEProblem(f,u0,(0.0,100.0))
```

However, this problem is supposed to conserve energy, and thus we define our manifold
to conserve the sum of squares:

```julia
function g(resid,u,p,t)
  resid[1] = u[2]^2 + u[1]^2 - 2
  resid[2] = 0
end
```

To build the callback, we just call

```julia
cb = ManifoldProjection(g)
```

Using this callback, the Runge-Kutta method `Vern7` conserves energy. Note that the
standard saving occurs after the step and before the callback, and thus we set
`save_everystep=false` to turn off all standard saving and let the callback
save after the projection is applied.

```julia
sol = solve(prob,Vern7(),save_everystep=false,callback=cb)
@test sol[end][1]^2 + sol[end][2]^2 ≈ 2
```

![manifold_projection](../assets/manifold_projection.png)

#### Saveat Warning

Note that the `ManifoldProjection` callback modifies the endpoints of the integration intervals
and thus breaks assumptions of internal interpolations. Because of this, the values for given by
saveat will not be order-matching. However, the interpolation error can be proportional to the
change by the projection, so if the projection is making small changes then one is still safe.
However, if there are large changes from each projection, you should consider only saving at
stopping/projection times. To do this, set `tstops` to the same values as `saveat`. There is a
performance hit by doing so because now the integrator is forced to stop at every saving point,
but this is guerenteed to match the order of the integrator even with the ManifoldProjection.

## AutoAbstol

Many problem solving environments [such as MATLAB](https://www.mathworks.com/help/simulink/gui/absolute-tolerance.html)
provide a way to automatically adapt the absolute tolerance to the problem. This
helps the solvers automatically "learn" what appropriate limits are. Via the
callback interface, DiffEqCallbacks.jl implements a callback `AutoAbstol` which
has the same behavior as the MATLAB implementation, that is the absolute tolerance
starts and at each iteration it is set to the maximum value that the state has thus
far reached times the relative tolerance. If `init_curmax` is zero, then the initial
value is determined by the `abstol` of the solver. Otherwise this is the initial
value for the current maximum `abstol`.

To generate the callback, use the constructor:

```julia
AutoAbstol(save=true;init_curmax=0.0)
```

## PositiveDomain

Especially in biology and other natural sciences, a desired property of
dynamical systems is the positive invariance of the positive cone, i.e.
non-negativity of variables at time $t_0$ ensures their non-negativity at times
$t \geq t_0$ for which the solution is defined. However, even if a system
satisfies this property mathematically it can be difficult for ODE solvers to
ensure it numerically, as these [MATLAB examples](https://www.mathworks.com/help/matlab/math/nonnegative-ode-solution.html)
show.

In order to deal with this problem one can specify `isoutofdomain=(u,p,t) -> any(x
-> x < 0, u)` as additional [solver option](@ref solver_options),
which will reject any step that leads to non-negative values and reduce the next
time step. However, since this approach only rejects steps and hence
calculations might be repeated multiple times until a step is accepted, it can
be computationally expensive.

Another approach is taken by a `PositiveDomain` callback in
DiffEqCallbacks.jl, which is inspired by
[Shampine's et al. paper about non-negative ODE solutions](http://www.sciencedirect.com/science/article/pii/S0096300304009683).
It reduces the next step by a certain scale factor until the extrapolated value
at the next time point is non-negative with a certain tolerance. Extrapolations
are cheap to compute but might be inaccurate, so if a time step is changed it
is additionally reduced by a safety factor of 0.9. Since extrapolated values are
only non-negative up to a certain tolerance and in addition actual calculations
might lead to negative values, also any negative values at the current time point
are set to 0. Hence by this callback non-negative values at any time point are
ensured in a computationally cheap way, but the quality of the solution
depends on how accurately extrapolations approximate next time steps.

Please note that the system should be defined also outside the positive domain,
since even with these approaches negative variables might occur during the
calculations. Moreover, one should follow Shampine's et. al. advice and set the
derivative ``x'_i`` of a negative component ``x_i`` to ``\max \{0, f_i(x, t)\}``,
where ``t`` denotes the current time point with state vector ``x`` and ``f_i``
is the ``i``-th component of function ``f`` in an ODE system ``x' = f(x, t)``.

### Constructor

```julia
PositiveDomain(u=nothing; save=true, abstol=nothing, scalefactor=nothing)
```

- `u`: A prototype of the state vector of the integrator. A copy of it is saved and
  extrapolated values are written to it. If it is not specified
  every application of the callback allocates a new copy of the state vector.
- `save`: Whether to do the standard saving (applied after the callback).
- `abstol`: Tolerance up to which negative extrapolated values are accepted.
  Element-wise tolerances are allowed. If it is not specified every application
  of the callback uses the current absolute tolerances of the integrator.
- `scalefactor`: Factor by which an unaccepted time step is reduced. If it is not
  specified time steps are halved.

## GeneralDomain

A `GeneralDomain` callback in DiffEqCallbacks.jl generalizes the concept of
a `PositiveDomain` callback to arbitrary domains. Domains are specified by
in-place functions `g(u, resid)` or `g(t, u, resid)` that calculate residuals of a
state vector `u` at time `t` relative to that domain. As for `PositiveDomain`, steps
are accepted if residuals of the extrapolated values at the next time step are below
a certain tolerance. Moreover, this callback is automatically coupled with a
`ManifoldProjection` that keeps all calculated state vectors close to the desired
domain, but in contrast to a `PositiveDomain` callback the nonlinear solver in a
`ManifoldProjection` can not guarantee that all state vectors of the solution are
actually inside the domain. Thus a `PositiveDomain` callback should in general be
preferred.

### Constructor

```julia
function GeneralDomain(g, u=nothing; nlsolve=NLSOLVEJL_SETUP(), save=true,
                       abstol=nothing, scalefactor=nothing, autonomous=numargs(g)==2,
                       nlopts=Dict(:ftol => 10*eps()))
```

- `g`: The residual function for the domain. This is an inplace function of form
  `g(resid, u, p, t)` which writes to the residual the difference from
  the domain.
- `u`: A prototype of the state vector of the integrator and the residuals. Two
  copies of it are saved, and extrapolated values and residuals are written to them.
  If it is not specified every application of the callback allocates two new copies
  of the state vector.
- `nlsolve`: A nonlinear solver as defined [in the nlsolve format](@ref linear_nonlinear)
  which is passed to a `ManifoldProjection`.
- `save`: Whether to do the standard saving (applied after the callback).
- `abstol`: Tolerance up to which residuals are accepted. Element-wise tolerances
  are allowed. If it is not specified every application of the callback uses the
  current absolute tolerances of the integrator.
- `scalefactor`: Factor by which an unaccepted time step is reduced. If it is not
  specified time steps are halved.
- `autonomous`: Whether `g` is an autonomous function of the form `g(u, resid)`.
- `nlopts`: Optional arguments to nonlinear solver of a `ManifoldProjection` which
  can be any of the [NLsolve keywords](https://github.com/JuliaNLSolvers/NLsolve.jl#fine-tunings).
  The default value of `ftol = 10*eps()` ensures that convergence is only declared
  if the infinite norm of residuals is very small and hence the state vector is very
  close to the domain.

## Stepsize Limiters

In many cases there is a known maximal stepsize for which the computation is
stable and produces correct results. For example, in hyperbolic PDEs one normally
needs to ensure that the stepsize stays below some ``\Delta t_{FE}`` determined
by the CFL condition. For nonlinear hyperbolic PDEs this limit can be a function
`dtFE(u,p,t)` which changes throughout the computation. The stepsize limiter lets
you pass a function which will adaptively limit the stepsizes to match these
constraints.

### Constructor

```julia
StepsizeLimiter(dtFE;safety_factor=9//10,max_step=false,cached_dtcache=0.0)
```

- `dtFE`: The function for the maximal timestep, called as `dtFE(u,p,t)`
  using the previous values of `u`, `p`, and `t`.
- `safety_factor`: The factor below the true maximum that will be stepped to
  which defaults to `9//10`.
- `max_step`: Makes every step equal to `safety_factor*dtFE(u,p,t)` when the
  solver is set to `adaptive=false`.
- `cached_dtcache`: Should be set to match the type for time when not using
  Float64 values.

## FunctionCallingCallback

The function calling callback lets you define a function `func(u,t,integrator)`
which gets calls at the time points of interest. The constructor is:

```julia
  FunctionCallingCallback(func;
                 funcat=Vector{Float64}(),
                 func_everystep=isempty(funcat),
                 func_start = true,
                 tdir=1)
```

- `func(u, t, integrator)` is the function to be called.
- `funcat` values that the function is sure to be evaluated at.
- `func_everystep` whether to call the function after each integrator step.
- `func_start` whether the function is called the initial condition.
- `tdir` should be `sign(tspan[end]-tspan[1])`. It defaults to `1` and should
    be adapted if `tspan[1] > tspan[end]`.

## SavingCallback

The saving callback lets you define a function `save_func(u, t, integrator)` which
returns quantities of interest that shall be saved.

### Constructor

```julia
SavingCallback(save_func, saved_values::SavedValues;
               saveat=Vector{eltype(saved_values.t)}(),
               save_everystep=isempty(saveat),
               tdir=1)
```
- `save_func(u, t, integrator)` returns the quantities which shall be saved.
  Note that this should allocate the output (not as a view to `u`).
- `saved_values::SavedValues` is the types that `save_func` will return, i.e.
  `save_func(u, t, integrator)::savevalType`. It's specified via
  `SavedValues(typeof(t),savevalType)`, i.e. give the type for time and the
  type that `save_func` will output (or higher compatible type).
- `saveat` mimicks `saveat` in `solve` from `solve`.
- `save_everystep` mimicks `save_everystep` from `solve`.
- `save_start` mimicks `save_start` from `solve`.
- `tdir` should be `sign(tspan[end]-tspan[1])`. It defaults to `1` and should
  be adapted if `tspan[1] > tspan[end]`.

The outputted values are saved into `saved_values`. Time points are found
via `saved_values.t` and the values are `saved_values.saveval`.

### Example

In this example we will solve a matrix equation and at each step save a tuple
of values which contains the current trace and the norm of the matrix. We build
the `SavedValues` cache to use `Float64` for time and `Tuple{Float64,Float64}`
for the saved values, and then call the solver with the callback.

```julia
using DiffEqCallbacks, OrdinaryDiffEq, LinearAlgebra
prob = ODEProblem((du,u,p,t) -> du .= u, rand(4,4), (0.0,1.0))
saved_values = SavedValues(Float64, Tuple{Float64,Float64})
cb = SavingCallback((u,t,integrator)->(tr(u),norm(u)), saved_values)
sol = solve(prob, Tsit5(), callback=cb)

print(saved_values.saveval)
#=
Tuple{Float64,Float64}[(2.23186, 2.49102), (2.46675, 2.75318), (3.16138, 3.52847), (4.42011, 4.93337), (6.06683, 6.77129)]
=#
```

Note that the values are retrieved from the cache as `.saveval`, and the time points are found as
`.t`. If we want to control the saved times, we use `saveat` in the callback. The save controls like
`saveat` act analogously to how they act in the `solve` function.

```julia
saved_values = SavedValues(Float64, Tuple{Float64,Float64})
cb = SavingCallback((u,t,integrator)->(tr(u),norm(u)), saved_values, saveat=0.0:0.1:1.0)
sol = solve(prob, Tsit5(), callback=cb)
print(saved_values.saveval)
print(saved_values.t)

#=
Tuple{Float64,Float64}[(2.23186, 2.49102), (2.46659, 2.753), (2.726, 3.04254), (3.0127, 3.36253),
(3.32955, 3.71617), (3.67972, 4.107), (4.06672, 4.53893), (4.49442, 5.0163), (4.9671, 5.54387),
(5.48949, 6.12692), (6.06683, 6.77129)]
[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
=#
```

## PresetTimeCallback

`PresetTimeCallback` is a callback that adds callback `affect!` calls at preset
times. No playing around with `tstops` or anything is required: this callback
adds the triggers for you to make it automatic.

```julia
PresetTimeCallback(tstops,user_affect!;
                            initialize = DiffEqBase.INITIALIZE_DEFAULT,
                            kwargs...)
```

- `tstops`: the times for the `affect!` to trigger at.
- `user_affect!`: an `affect!(integrator)` function to use at the time points.

## IterativeCallback

`IterativeCallback` is a callback to be used to iteratively apply some effect.
For example, if given the first effect at `t₁`, you can define `t₂` to apply
the next effect.

A `IterativeCallback` is constructed as follows:

```julia
function IterativeCallback(time_choice, user_affect!,tType = Float64;
                           initial_affect = false, kwargs...)
```

where `time_choice(integrator)` determines the time of the next callback and
`user_affect!` is the effect applied to the integrator at the stopping points.
If `nothing` is returned for the time choice then the iterator ends. `initial_affect`
is whether to apply the affect at `t=0` which defaults to `false`. `kwargs` are
keyword arguments accepted by the `DiscreteCallback` constructor.

## PeriodicCallback

`PeriodicCallback` can be used when a function should be called periodically in
terms of integration time (as opposed to wall time), i.e. at `t = tspan[1]`,
`t = tspan[1] + Δt`, `t = tspan[1] + 2Δt`, and so on. This callback can, for
example, be used to model a digital controller for an analog system, running at
a fixed rate.

### Constructor

```julia
PeriodicCallback(f, Δt::Number; initial_affect = true, kwargs...)
```

where `f` is the function to be called periodically, `Δt` is the period,
`initial_affect` is whether to apply the affect at `t=0` which defaults to `true`,
and `kwargs` are keyword arguments accepted by the `DiscreteCallback` constructor
(see the [DiscreteCallback](@ref) section).

## TerminateSteadyState

`TerminateSteadyState` can be used to solve the problem for the steady-state
by running the solver until the derivatives of the problem converge to 0 or
`tspan[2]` is reached. This is an alternative approach to root finding (see
the [Steady State Solvers](@ref) section). The constructor of this callback is:

```julia
TerminateSteadyState(abstol = 1e-8, reltol = 1e-6, test = allDerivPass)
```

where `abstol` and `reltol` are the absolute and relative tolerance, respectively.
These tolerances may be specified as scalars or as arrays of the same length
as the states of the problem. `test` represents the function that evaluates the
condition for termination. The default condition is that all derivatives should
become smaller than `abstol` and the states times `reltol`. The user
can pass any other function to implement a different termination condition. Such
function should take four arguments: `integrator` (see [Integrator Interface](@ref)
for details), `abstol` and `reltol`.
