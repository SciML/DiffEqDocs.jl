# [Integrator Interface](@id integrator)

The integrator interface gives one the ability to interactively step through
the numerical solving of a differential equation. Through this interface,
one can easily monitor results, modify the problem during a run, and dynamically
continue solving as one sees fit.

## Initialization and Stepping

To initialize an integrator, use the syntax:

```julia
integrator = init(prob, alg; kwargs...)
```

The keyword args which are accepted are the same as the [solver options](@ref solver_options)
used by `solve` and the returned value is an `integrator` which satisfies
`typeof(integrator)<:DEIntegrator`. One can manually choose to step via the `step!` command:

```julia
step!(integrator)
```

which will take one successful step. Additionally:

```julia
step!(integrator, dt, false)
```

passing a `dt` will make the integrator continue to step until at least `integrator.t+dt`, and
passing `true` as the third argument will add a `tstop` to force it to step up to `integrator.t+dt`,
exactly.

To check whether the integration step was successful, you can
call `check_error(integrator)` which returns one of the
[return codes](@ref retcodes).

This type also implements an iterator interface, so one can step `n` times
(or to the last `tstop`) using the `take` iterator:

```julia
for i in take(integrator, n)
end
```

One can loop to the end by using `solve!(integrator)` or using the iterator interface:

```julia
for i in integrator
end
```

In addition, some helper iterators are provided to help monitor the solution. For
example, the `tuples` iterator lets you view the values:

```julia
for (u, t) in tuples(integrator)
    @show u, t
end
```

and the `intervals` iterator lets you view the full interval:

```julia
for (uprev, tprev, u, t) in intervals(integrator)
    @show tprev, t
end
```

Additionally, you can make the iterator return specific time points via the
`TimeChoiceIterator`:

```julia
ts = range(0, stop = 1, length = 11)
for (u, t) in TimeChoiceIterator(integrator, ts)
    @show u, t
end
```

Lastly, one can dynamically control the “endpoint”. The initialization simply makes
`prob.tspan[2]` the last value of `tstop`, and many of the iterators are made to stop
at the final `tstop` value. However, `step!` will always take a step, and one
can dynamically add new values of `tstops` by modifying the variable in the
options field: `add_tstop!(integrator,new_t)`.

Finally, to solve to the last `tstop`, call `solve!(integrator)`. Doing `init`
and then `solve!` is equivalent to `solve`.

```@docs
SciMLBase.step!
SciMLBase.check_error
SciMLBase.check_error!
```

## Handling Integrators

The `integrator<:DEIntegrator` type holds all the information for the intermediate solution
of the differential equation. Useful fields are:

  - `t` - time of the proposed step
  - `u` - value at the proposed step
  - `p` - user-provided data
  - `opts` - common solver options
  - `alg` - the algorithm associated with the solution
  - `f` - the function being solved
  - `sol` - the current state of the solution
  - `tprev` - the last timepoint
  - `uprev` - the value at the last timepoint
  - `tdir` - the sign for the direction of time

The function `f` is usually a wrapper of the function provided when creating the
specific problem. For example, when solving an `ODEProblem`, `f` will be an
`ODEFunction`. To access the right-hand side function provided by the user when
creating the `ODEProblem`, please use `SciMLBase.unwrapped_f(integrator.f.f)`.

The `p` is the (parameter) data which is provided by the user as a keyword arg in
`init`. `opts` holds all the common solver options, and can be mutated to
change the solver characteristics. For example, to modify the absolute tolerance
for the future timesteps, one can do:

```julia
integrator.opts.abstol = 1e-9
```

The `sol` field holds the current solution. This current solution includes the
interpolation function if available, and thus `integrator.sol(t)` lets one
interpolate efficiently over the whole current solution. Additionally,
a “current interval interpolation function” is provided on the `integrator` type
via `integrator(t,deriv::Type=Val{0};idxs=nothing,continuity=:left)`.
This uses only the solver information from the interval
`[tprev,t]` to compute the interpolation, and is allowed to extrapolate beyond
that interval.

### Note about mutating

Be cautious: one should not directly mutate the `t` and `u` fields of the integrator.
Doing so will destroy the accuracy of the interpolator and can harm certain algorithms.
Instead, if one wants to introduce discontinuous changes, one should use the
[callbacks](@ref callbacks). Modifications within a callback
`affect!` surrounded by saves provides an error-free handling of the discontinuity.

As low-level alternative to the callbacks, one can use `set_t!`, `set_u!` and
`set_ut!` to mutate integrator states.  Note that certain integrators may not
have efficient ways to modify `u` and `t`.  In such case, `set_*!` are as
inefficient as `reinit!`.

```@docs
SciMLBase.set_t!
SciMLBase.set_u!
SciMLBase.set_ut!
```

### Integrator vs Solution

The integrator and the solution have very different actions because they have
very different meanings. The `typeof(sol) <: DESolution` type is a type with
history: it stores
all the (requested) timepoints and interpolates/acts using the values closest
in time. On the other hand, the `typeof(integrator)<:DEIntegrator` type is a
local object. It only knows the times of the interval it currently spans,
the current caches and values, and the current state of the solver
(the current options, tolerances, etc.). These serve very different purposes:

  - The `integrator`'s interpolation can extrapolate, both forward and backward
    in time. This is used to estimate events and is internally used for predictions.
  - The `integrator` is fully mutable upon iteration. This means that every time
    an iterator affect is used, it will take timesteps from the current time. This
    means that `first(integrator)!=first(integrator)` since the `integrator` will
    step once to evaluate the left and then step once more (not backtracking).
    This allows the iterator to keep dynamically stepping, though one should note
    that it may violate some immutability assumptions commonly made about iterators.

If one wants the solution object, then one can find it in `integrator.sol`.

## Function Interface

In addition to the type interface, a function interface is provided which allows
for safe modifications of the integrator type, and allows for uniform usage
throughout the ecosystem (for packages/algorithms which implement the functions).
The following functions make up the interface:

### Saving Controls

```@docs
savevalues!
```

### Caches

```@docs
get_tmp_cache
full_cache
```

### [Stepping Controls](@id stepping_controls)

```@docs
u_modified!
get_proposed_dt
set_proposed_dt!
terminate!
change_t_via_interpolation!
add_tstop!
has_tstop
first_tstop
pop_tstop!
add_saveat!
```

### Resizing

```@docs
resize!
deleteat!
addat!
resize_non_user_cache!
deleteat_non_user_cache!
addat_non_user_cache!
```

### Reinitialization

```@docs
reinit!
auto_dt_reset!
```

### Misc

```@docs
get_du
get_du!
```

!!! warning
    
    Note that not all these functions will be implemented for every algorithm.
    Some have hard limitations. For example, Sundials.jl cannot resize problems.
    When a function is not limited, an error will be thrown.

## Additional Options

The following options can additionally be specified in `init` (or be mutated in
the `opts`) for further control of the integrator:

  - `advance_to_tstop`: This makes `step!` continue to the next value in `tstop`.
  - `stop_at_next_tstop`: This forces the iterators to stop at the next value of `tstop`.

For example, if one wants to iterate but only stop at specific values, one can
choose:

```julia
integrator = DE.init(
    prob, DE.Tsit5(); dt = 1 // 2^(4), tstops = [0.5], advance_to_tstop = true)
for (u, t) in tuples(integrator)
    @test t ∈ [0.5, 1.0]
end
```

which will only enter the loop body at the values in `tstops` (here, `prob.tspan[2]==1.0`
and thus there are two values of `tstops` which are hit). Additionally, one can
`solve!` only to `0.5` via:

```julia
integrator = DE.init(prob, DE.Tsit5(); dt = 1 // 2^(4), tstops = [0.5])
integrator.opts.stop_at_next_tstop = true
solve!(integrator)
```

## Plot Recipe

Like the `DESolution` type, a plot recipe is provided for the `DEIntegrator` type.
Since the `DEIntegrator` type is a local state type on the current interval,
`plot(integrator)` returns the solution on the current interval. The same
options for the plot recipe are provided as for `sol`, meaning one can choose
variables via the `idxs` keyword argument, or change the `plotdensity` / turn
on/off `denseplot`.

Additionally, since the `integrator` is an iterator, this can be used in the
Plots.jl `animate` command to iteratively build an animation of the solution
while solving the differential equation.

For an example of manually chaining together the iterator interface and plotting,
one should try the following:

```julia
import DifferentialEquations as DE, DiffEqProblemLibrary, Plots

# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
prob = DE.ODEProblem((u, p, t) -> 1.01u, 0.5, (0.0, 1.0))

import Plots
integrator = DE.init(prob, DE.Tsit5(); dt = 1 // 2^(4), tstops = [0.5])
pyplot(show = true)
Plots.plot(integrator)
for i in integrator
    display(Plots.plot!(integrator, idxs = (0, 1), legend = false))
end
DE.step!(integrator);
Plots.plot!(integrator, idxs = (0, 1), legend = false);
savefig("iteratorplot.png")
```

![Iterator Plot](../assets/iteratorplot.png)
