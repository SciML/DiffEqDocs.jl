# Integrator Interface

The integrator interface gives one the ability to interactively step through
the numerical solving of a differential equation. Through this interface,
one can easily monitor results, modify the problem during a run, and dynamically
continue solving as one sees fit.

## Initialization and Stepping

To initialize an integrator, use the syntax:

```julia
integrator = init(prob,alg;kwargs...)
```

The keyword args which are accepted are the same [Common Solver Options](@ref)
used by `solve`. The type which is returned is the integrator. One can manually
choose to step via the `step!` command:

```julia
step!(integrator)
```

which will take one successful step. This type also implements an integrator interface,
so one can step `n` times (or to the last `tstop`) using the `take` iterator:

```julia
for i in take(integrator,n) end
```

One can loop to the end by using `solve!(integrator)` or using the integrator interface:

```julia
for i in integrator end
```

In addition, some helper iterators are provided to help monitor the solution. For
example, the `tuples` iterator lets you view the values:

```julia
for (t,u) in tuples(integrator)
  @show t,u
end
```

and the `intervals` iterator lets you view the full interval:

```julia
for (tprev,uprev,t,u) in intervals(integrator)
  @show tprev,t
end
```

Lastly, one can dynamically control the "endpoint". The initialization simply makes
`prob.tspan[2]` the last value of `tstop`, and many of the iterators are made to stop
at the final `tstop` value. However, `step!` will always take a step, and one
can dynamically add new values of `tstops` by modifiying the variable in the
options field: `add_tstop!(integrator,new_t)`.

## Handing Integrators

The `integrator` type holds all of the information for the intermediate solution
of the differential equation. Useful fields are:

* `t` - time of the proposed step
* `u` - value at the proposed step
* `userdata` - user-provided data type
* `opts` - common solver options
* `alg` - the algorithm associated with the solution
* `f` - the function being solved
* `sol` - the current state of the solution
* `tprev` - the last timepoint
* `uprev` - the value at the last timepoint

The `userdata` is the type which is provided by the user as a keyword arg in
`init`. `opts` holds all of the common solver options, and can be mutated to
change the solver characteristics. For example, to modify the absolute tolerance
for the future timesteps, one can do:

```julia
integrator.opts.abstol = 1e-9
```

The `sol` field holds the current solution. This current solution includes the
interpolation function if available, and thus `integrator.sol(t)` lets one
interpolate efficiently over the whole current solution. Additionally, a
a "current interval interpolation function" is provided on the `integrator` type
via `integrator(t)`. This uses only the solver information from the interval
`[tprev,t]` to compute the interpolation, and is allowed to extrapolate beyond
that interval.

### Note about mutating

Be cautious: one should not directly mutate the `t` and `u` fields of the integrator.
Doing so will destroy the accuracy of the interpolator and can harm certain algorithms.
Instead if one wants to introduce discontinuous changes, one should use the
[Event Handling and Callback Functions](@ref). Modifications within a callback
`affect!` surrounded by saves provides an error-free handling of the discontinuity.

### Integrator vs Solution

The integrator and the solution have very different actions because they have
very different meanings. The `Solution` type is a type with history: it stores
all of the (requested) timepoints and interpolates/acts using the values closest
in time. On the other hand, the `Integrator` type is a local object. It only knows
the times of the interval it currently spans, the current caches and values,
and the current state of the solver (the current options, tolerances, etc.).
These serve very different purposes:

* The `integrator`'s interpolation can extrapolate, both forward and backward in
  in time. This is used to estimate events and is internally used for predictions.
* The `integrator` is fully mutable upon iteration. This means that every time
  an iterator affect is used, it will take timesteps from the current time. This
  means that `first(integrator)!=first(integrator)` since the `integrator` will
  step once to evaluate the left and then step once more (not backtracking).
  This allows the iterator to keep dynamically stepping, though one should note
  that it may violate some immutablity assumptions commonly made about iterators.

If one wants the solution object, then one can find it in `integrator.sol`.

## Function Interface

In addition to the type interface, a function interface is provided which allows
for safe modifications of the integrator type, and allows for uniform usage
throughout the ecosystem (for packages/algorithms which implement the functions).
The following functions make up the interface:

### Saving Controls

* `savevalues!(integrator)`: Adds the current state to the `sol`.

### Caches

* `get_tmp_cache(integrator)`: Returns a tuple of internal cache vectors which are
  safe to use as temporary arrays. This should be used for integrator interface
  and callbacks which need arrays to write into in order to be non-allocating.
  The length of the tuple is dependent on the method.
* `user_cache(integrator)`: Returns an iterator over the user-facing cache arrays.
* `u_cache(integrator)`:  Returns an iterator over the cache arrays for `u` in the method.
  This can be used to change internal values as needed.
* `du_cache(integrator)`:  Returns an iterator over the cache arrays for rate quantities the method.
  This can be used to change internal values as needed.
* `full_cache(integrator)`:  Returns an iterator over the cache arrays of the method.
  This can be used to change internal values as needed.

### Stepping Controls

* `u_modified!(integrator,bool)`: Bool which states whether a change to `u` occurred,
  allowing the solver to handle the discontinuity. By default, this is assumed
  to be true if a callback is used. This will result in the re-calculation of
  the derivative at `t+dt`, which is not necessary if the algorithm is FSAL
  and `u` does not experience a discontinuous change at the end of the interval.
  Thus if `u` is unmodified in a callback, a single call to the derivative calculation
  can be eliminated by `u_modified!(integrator,false)`.
* `get_proposed_dt(integrator)`:  Gets the proposed `dt` for the
  next timestep.
* `set_proposed_dt!(integrator,dt)`:  Sets the proposed `dt` for the
  next timestep.
* `set_proposed_dt!(integrator,integrator2)`:  Sets the timestepping of
  `integrator` to match that of `integrator2`. Note that due to PI control and
  step acceleration this is more than matching the factors in most cases.
* `proposed_dt(integrator)`: Returns the `dt` of the proposed step.
* `terminate!(integrator)`: Terminates the integrator by emptying `tstops`. This
  can be used in events and callbacks to immediately end the solution process.
* `change_t_via_interpolation(integrator,t,modify_save_endpoint=Val{false})`: This
  option lets one modify the current `t` and changes all of the corresponding
  values using the local interpolation. If the current solution has already
  been saved, one can provide the optional value `modify_save_endpoint` to also
  modify the endpoint of `sol` in the same manner.
* `add_tstop!(integrator,t)`: Adds a `tstop` at time `t`.
* `add_saveat!(integrator,t)`: Adds a `saveat` time point at `t`.

### Resizing

* `resize!(integrator,k)`: Resizes the DE to a size `k`. This chops off the end
  of the array, or adds blank values at the end, depending on whether `k>length(integrator.u)`.
* `resize_non_user_cache!(integrator,k)`: Resizes the non-user facing caches to be
  compatible with a DE of size `k`. This includes resizing Jacobian caches. Note
  that in many cases, `resize!` simple resizes `user_cache` variables and then
  calls this function. This finer control is required for some `AbstractArray`
  operations.
* `deleteat_non_user_cache!(integrator,idxs)`: `deleteat!`s the non-user facing caches
  at indices `idxs`. This includes resizing Jacobian caches. Note
  that in many cases, `deleteat!` simple `deleteat!`s `user_cache` variables and then
  calls this function. This finer control is required for some `AbstractArray`
  operations.
* `addat_non_user_cache!(integrator,idxs)`: `addat!`s the non-user facing caches
  at indices `idxs`. This includes resizing Jacobian caches. Note
  that in many cases, `addat!` simple `addat!`s `user_cache` variables and then
  calls this function. This finer control is required for some `AbstractArray`
  operations.
* `deleteat!(integrator,idxs)`: Shrinks the ODE by deleting the `idxs` components.
* `addat!(integrator,idxs)`: Grows the ODE by adding the `idxs` components.
  Must be contiguous indices.

### Reinit

The reinit function lets you restart the integration at a new value. The full
function is of the form:

```julia
reinit!(integrator::ODEIntegrator,u0 = integrator.sol.prob.u0;
  t0 = integrator.sol.prob.tspan[1], tf = integrator.sol.prob.tspan[2],
  erase_sol = true,
  tstops = integrator.opts.tstops_cache,
  saveat = integrator.opts.saveat_cache,
  d_discontinuities = integrator.opts.d_discontinuities_cache,
  reset_dt = (integrator.dtcache == zero(integrator.dt)) && integrator.opts.adaptive,
  reinit_callbacks = true, initialize_save = true,
  reinit_cache = true)
```

`u0` is the value to start at. The starting time point and end point can be changed
via `t0` and `tf`. `erase_sol` allows one to start with no other values in the
solution, or keep the previous solution. `tstops`, `d_discontinuities`, and
`saveat` are reset as well, but can be ignored. `reset_dt` is a boolean for
whether to reset the current value of `dt` using the automatic `dt` determination
algorithm. `reinit_callbacks` is whether to run the callback initializations
again (and `initialize_save` is for that). `reinit_cache` is whether to re-run
the cache initialization function (i.e. resetting FSAL, not allocating vectors)
which should usually be true for correctness.

Additionally, once can access `auto_dt_reset!(integrator::ODEIntegrator)` which
will run the auto `dt` initialization algorithm. 

### Misc

* `get_du(integrator)`: Returns the derivative at `t`.
* `get_du!(out,integrator)`: Write the current derivative at `t` into `out`.

#### Note

Note that not all of these functions will be implemented for every algorithm.
Some have hard limitations. For example, Sundials.jl cannot resize problems.
When a function is not limited, an error will be thrown.

## Additional Options

The following options can additionally be specified in `init` (or be mutated in
the `opts`) for further control of the integrator:

* `advance_to_tstop`: This makes `step!` continue to the next value in `tstop`.
* `stop_at_next_tstop`: This forces the iterators to stop at the next value of `tstop`.

For example, if one wants to iterate but only stop at specific values, one can
choose:

```julia
integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5],advance_to_tstop=true)
for (t,u) in tuples(integrator)
  @test t ∈ [0.5,1.0]
end
```

which will only enter the loop body at the values in `tstops` (here, `prob.tspan[2]==1.0`
and thus there are two values of `tstops` which are hit). Addtionally, one can
`solve!` only to `0.5` via:

```julia
integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])
integrator.opts.stop_at_next_tstop = true
solve!(integrator)
```

## Plot Recipe

Like the `Solution` type, a plot recipe is provided for the `Integrator` type.
Since the `Integrator` type is a local state type on the current interval,
`plot(integrator)` returns the solution on the current interval. The same
options for the plot recipe are provided as for `sol`, meaning one can choose
variables via the `vars` keyword argument, or change the `plotdensity` / turn
on/off `denseplot`.

Additionally, since the `integrator` is an integrator, this can be used in the
Plots.jl `animate` command to iteratively build an animation of the solution
while solving the differential equation.

For an example of manually chaining together the iterator interface and plotting,
one should try the following:

```julia
using DifferentialEquations, DiffEqProblemLibrary, Plots

# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
prob = ODEProblem((t,u)->1.01u,0.5,(0.0,1.0))

using Plots
integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])
pyplot(show=true)
plot(integrator)
for i in integrator
  display(plot!(integrator,vars=(0,1),legend=false))
end
step!(integrator); plot!(integrator,vars=(0,1),legend=false)
savefig("iteratorplot.png")
```

![Iterator Plot](../assets/iteratorplot.png)
