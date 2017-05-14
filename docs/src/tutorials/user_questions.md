User case walkthrough

A user comes into chat and asks:

"I know SSRK methods are already available, so I would like to reuse them.
 I have been researching the diffEqs developer guide, however I am relatively new on Julia.
 I have developed similar methods on Matlab, but I am interested on using Julia for my work, so it seem a good point to start.
 I would appreciate any help.
 Maybe someone can share some links of related code that I can review or give me general advice and tips.
 Thanks in advance."

A good place to start is to know that the overarching logic is handled by two things:
This [initialization](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/blob/master/src/solve.jl#L11)

```
function init{algType<:OrdinaryDiffEqAlgorithm,recompile_flag}(
  prob::AbstractODEProblem,
  alg::algType,timeseries_init=typeof(prob.u0)[],ts_init=eltype(prob.tspan)[],ks_init=[],
  recompile::Type{Val{recompile_flag}}=Val{true};
  timeseries_steps = 1,
  saveat = eltype(prob.tspan)[],tstops = eltype(prob.tspan)[],d_discontinuities= eltype(prob.tspan)[],
  save_idxs = nothing,
  save_everystep = isempty(saveat),
  save_timeseries = nothing,save_start = true,
  dense = save_everystep && !(typeof(alg) <: Discrete),
  calck = (!isempty(setdiff(saveat,tstops)) || dense),
  dt = typeof(alg) <: Discrete && isempty(tstops) ? eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
  adaptive = isadaptive(alg),
  gamma=9//10,
  abstol=nothing,
  reltol=nothing,
  qmax=qmax_default(alg),qmin=qmin_default(alg),
  qoldinit=1//10^4, fullnormalize=true,
  beta2=beta2_default(alg),
  beta1=beta1_default(alg,beta2),
  maxiters = 1000000,
  dtmax=eltype(prob.tspan)((prob.tspan[end]-prob.tspan[1])),
  dtmin=eltype(prob.tspan) <: AbstractFloat ? eltype(prob.tspan)(10)*eps(eltype(prob.tspan)) : eltype(prob.tspan)(1//10^(10)),
  internalnorm = ODE_DEFAULT_NORM,
  isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
  unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
  verbose = true, force_dtmin = false,
  timeseries_errors = true, dense_errors=false,
  advance_to_tstop = false,stop_at_next_tstop=false,
  progress=false,progress_steps=1000,progress_name="ODE",
  progress_message = ODE_DEFAULT_PROG_MESSAGE,
  userdata=nothing,callback=nothing,
  allow_extrapolation = alg_extrapolates(alg),
  initialize_integrator=true,kwargs...)
  ```

  and this [solve!](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/blob/master/src/solve.jl#L303)

  ```

function solve!(integrator::ODEIntegrator)
  @inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      @ode_exit_conditions
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
      if isempty(integrator.opts.tstops)
        break
      end
    end
    handle_tstop!(integrator)
  end

  ```

  As such, there is no need to focus on the main details and instead define those
  two things.

  Now, an algorithm is defined by its cache and its `perform_step`.

  Let's look at Euler as an example here. The cache holds the arrays that [are updated:](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/blob/master/src/caches.jl#L37)
  ```
  immutable EulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
   u::uType
   uprev::uType
   tmp::uType
   k::rateType
   fsalfirst::rateType
 end

```

and performs the Euler steps, either on a [numbers:](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/blob/master/src/integrators/fixed_timestep_integrators.jl#L42)
```
@inline function perform_step!(integrator,cache::EulerConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  k = integrator.fsalfirst
  u = muladd(dt,k,uprev)
  k = f(t+dt,u) # For the interpolation, needs k at the updated point
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  @pack integrator = t,dt,u
end
```

or on [arrays:](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/blob/master/src/integrators/fixed_timestep_integrators.jl#L64)
```
@inline function perform_step!(integrator,cache::EulerCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  for i in uidx
    u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
  end
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  @pack integrator = t,dt,u,k
end
```
You may have to initialize before the first step, but that's just filling in the cache
arrays and setting up [pointers:](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/blob/master/src/integrators/fixed_timestep_integrators.jl#L53)
```
@inline function initialize!(integrator,cache::EulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # For the interpolation, needs k at the updated point
end
```

But that's it: everything else is generic code or other algorithms, though we
should take not that Euler is writtern peculiarly:

```
@inline function perform_step!(integrator,cache::EulerCache,f=integrator.f)
  @unpack t,dt,uprev,u,k = integrator
  uidx = eachindex(integrator.uprev)
  for i in uidx
    u[i] = muladd(dt,integrator.fsalfirst[i],uprev[i])
  end
  f(t+dt,u,integrator.fsallast) # For the interpolation, needs k at the updated point
  @pack integrator = t,dt,u,k
end
```

because this way allows the derivative to be calculated at both endpoints efficiently so that way an interpolation exists.
But the algorith is just instead of
$u_{n+1} = u_n + Δtf(t,u_n)$
you do
$f_n = f(t,u_n)$
to initialize, and then each step is
$u_{n+1} = u_n + Δtf_n, f_{n+1} = f(t,u_{n+1})$
(since ($f_n, f_{n+1},u_n, u_{n+1}$) gives you a 3rd order interpolation of the interval.)
Although those are extra details to be handled later - your first information can
ignore all of that First Same As Last (FSAL) stuff and just calculate `u` from `uprev`.

"Thanks! . On the step function I would like to use Runge Kutta for time... I can call other integrator, or is better just copy the code inside the step function?"

 It might be easier to copy over some code for now, but a way to handle this would be to have a Runge-Kutta method's cache as part of your method's cache, and yes call `perform_step` on that inside of your `perform_step`.

 Go for simple implementations first, and build up the Julia type-system knowledge for full exploitation of the method later.

 "For the additional parameters required by the method (like CFL condition or the Jacobian of f(x)), I have to setup a new problem type...
or I just include these as options options on the initialize function..."

Jacobians are handled through [dispatches:](http://docs.juliadiffeq.org/latest/features/performance_overloads.html#Declaring-Explicit-Jacobians-1), but they will soon be handled by @miguelraz's GSoc work. 

The Rosenbrock code is a bit too contrived for a first example of autodifferentiation, but you may look at it if needed.

For the CFL condition, you can require extra parameters provided to the algorithm type itself.

An algorithm without default can happen, as in `sol = solve(prob,MyAlgorithm())`, but with CFL defined it could `sol = solve(prob, MyAlgorithm(CFL=1.0))`.
