# Common Solver Options

The DifferentialEquations.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. Not all of the interface is provided by every algorithm.
For more detailed information on the defaults and the available options
for specific algorithms / packages, see the manual pages for the solvers of specific
problems.

## Default Algorithm Hinting

To help choose the default algorithm, the keyword argument `alg_hints` is provided.
`alg_hints` is a `Vector{Symbol}` which describe the problem at a high level
to the solver. The options are:

* `:nonstiff` - Denotes the equation as nonstiff.
* `:stiff` - Denotes the equation as stiff.

Currently unused options include:

* `:interpolant` - Denotes that a high-precision interpolation is important.
* `:memorybound` - Denotes that the solver will be memory bound.

This functionality is derived via the benchmarks in [DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl)
and is under active development.

## Output Control

These arguments control the output behavior the solvers. It defaults to maximum
output to give the best interactive user experience, but can be reduced all the
way to only saving the solution at the final timepoint. All of these options
can be mixed and matched. For example, the combination

```julia
saveat=[0:1/4:1],save_timeseries=false,dense=false
```

will only save the solution at the timepoints `0:1/4:1` and no other locations.
For more examples for controlling the output behavior, see the
[Output Specification manual page](../man/output_specification.html).

* `dense`: Denotes whether to save the extra pieces for dense (continuous)
  output. Default is true for algorithms which have the ability to produce dense output.
* `saveat`: Denotes extra times to save the solution at during the solving
  phase. The solver will save at each of the timepoints in this array in the
  most efficient manner. Note that this can be used even if `dense=false`.
  For methods where interpolation is not possible, this may be equivalent to
  `tstops`. Default is `[]`.
* `tstops`: Denotes extra times that the timestepping algorithm must step to.
  This should be used to help the solver deal with discontinuities and
  singularities, since stepping exactly at the time of the discontinuity will
  improve accuracy. If a method cannot change timesteps (fixed timestep multistep methods),
  then `tstops` will use an interpolation, matching the behavior of `saveat`.
  If a method cannot change timesteps and also cannot interpolation,
  then `tstops` must be a multiple of `dt` or elese an error will be thrown. Default is `[]`.
* `d_discontinuities:` Denotes locations of disccontinuities in low order derivatives.
  This will force FSAL algorithms which assume derivative continuity to re-evaluate
  the derivatives at point of discontinuity. The default is `[]`.
* `calck`: Turns on and off the internal ability for intermediate    
  interpolations. This defaults to `dense || !isempty(saveat) || `
  "no custom callback is given". This can be used to turn off interpolations
  (to save memory) if one isn't using interpolations when a custom callback is
  used. Another case where this may be used is to turn on interpolations for
  usage in the integrator interface even when interpolations are used no
   where else. Note that this is only required if the algorithm doesn't have
  a free or lazy interpolation (`DP8()`).
* `save_timeseries`: Saves the result at every timeseries_steps iteration.    
  Default is true.
* `timeseries_steps`: Denotes how many steps between saving a value for the
  timeseries. Defaults to 1.

## Stepsize Control

These arguments control the timestepping routines.

* `adaptive`: Turns on adaptive timestepping for appropriate methods. Default
  is true.
* `abstol`: Absolute tolerance in adaptive timestepping. Defaults to 1e-3.
* `reltol`: Relative tolerance in adaptive timestepping. Defaults to 1e-6.
* `dt`: Sets the initial stepsize. This is also the stepsize for fixed
  timestep methods. Defaults to an automatic choice.
* `internalnorm`: The norm function `internalnorm(u)` which error estimates
  are calculated.
  Defaults are package-dependent.
* `gamma`: The risk-factor γ in the q equation for adaptive timestepping.
  Default is algorithm dependent.
* `dtmax`: Maximum dt for adaptive timestepping. Defaults are
  package-dependent.
* `dtmin`: Minimum dt for adaptive timestepping. Defaults are
  package-dependent.
* `beta1`: The Lund stabilization α parameter. Defaults are
  algorithm-dependent.
* `beta2`: The Lund stabilization β parameter. Defaults are
  algorithm-dependent.
* `qmax`: Defines the maximum value possible for the adaptive q. Defaults are
  algorithm-dependent.
* `qmin`: Defines the maximum value possible for the adaptive q. Defaults are
  algorithm-dependent.
* `qoldinit`: The initial `qold` in stabilization stepping. Defaults are
  algorithm-dependent.

### Fixed Stepsize Usage

Note that if a method does not have adaptivity, the following rules apply:

* If `dt` is set, then the algorithm will step with size `dt` each iteration.
* If `tstops` and `dt` are both set, then the algorithm will step with either a
  size `dt`, or use a smaller step to hit the `tstops` point.
* If `tstops` is set without `dt`, then the algorithm will step directly to
  each value in `tstops`
* If neither `dt` nor `tstops` are set, the solver will throw an error.

## Miscellaneous

* `maxiters`: Maximum number of iterations before stopping. Defaults to 1e5.
* `callback`: Specifies a callback. Defaults to a callback function which
  performs the saving routine. For more information, see the
  [Event Handling and Callback Functions manual page](https://juliadiffeq.github.io/DiffEqDocs.jl/latest/man/callback_functions.html).
* `isoutofdomain`: Specifies a function `isoutofdomain(t,u)` where, when it
  returns false, it will reject the timestep. Defaults to always false.

## Progress Monitoring

These arguments control the usage of the progressbar in the Juno IDE.

* `progress`: Turns on/off the Juno progressbar. Default is false.
* `progress_steps`: Numbers of steps between updates of the progress bar.
  Default is 1000.
* `progress_name`: Controls the name of the progressbar. Default is the name
  of the problem type.
* `progress_message`: Controls the message with the progressbar. Defaults to
  showing `dt`, `t`, the maximum of `u`.

## User Data

* `userdata`: This is a user-chosen type which will show up in the `integrator`
  type, allowing the user to have a cache for callbacks, event handling, and
  other various activities.

## Error Calculations

If you are using the test problems (ex: `ODETestProblem`), then the following
options control the errors which are calculated:

* `timeseries_errors`: Turns on and off the calculation of errors at the steps
  which were taken, such as the `l2` error. Default is true.
* `dense_errors`: Turns on and off the calculation of errors at the steps which
  require dense output and calculate the error at 100 evenly-spaced points
  throughout `tspan`. An example is the `L2` error. Default is false.
