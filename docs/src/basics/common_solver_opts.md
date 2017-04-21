# Common Solver Options

The DifferentialEquations.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. Not all of the interface is provided by every algorithm.
For more detailed information on the defaults and the available options
for specific algorithms / packages, see the manual pages for the solvers of specific
problems. To see whether a specific package is compaible with the use of a
given option, see the [compatibility chart](compatibility_chart.html)

## Default Algorithm Hinting

To help choose the default algorithm, the keyword argument `alg_hints` is provided to `solve`.
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

These arguments control the output behavior of the solvers. It defaults to maximum
output to give the best interactive user experience, but can be reduced all the
way to only saving the solution at the final timepoint. All of these options
can be mixed and matched. For example, the combination:

```julia
sol = solve(prob; saveat=[0.2, 0.5], dense = true)
```

will only save the solution (`sol.u`) at the timepoints `tspan[1], 0.2, 0.5, tspan[end]`. 
It will also enable dense output to the `sol` object, enabling you to do something like `sol(0.345)` which interpolates
the solution to the time equal to 0.345.
For more examples for controlling the output behavior, see the
[Output Specification manual page](../features/output_specification.html).

* `dense`: Denotes whether to save the extra pieces required for dense (continuous)
  output. Default is true for algorithms which have the ability to produce dense output.
* `saveat`: Denotes specific times to save the solution at, during the solving
  phase. The solver will save at each of the timepoints in this array in the
  most efficient manner (always including the points of `tspan`). 
  Note that this can be used even if `dense=false`. In fact, if only `saveat` is given, then the 
  arguments `save_everystep` and `dense` are becoming `false` by default and must be explicitly given as `true` if desired.
  If `saveat` is given a number, then it will
  automatically expand to `tspan[1]:saveat:tspan[2]`. 
  For methods where interpolation is not possible, `saveat` may be equivalent to
  `tstops`. Default is `[]`. 
* `save_idxs`: Denotes the indices for the components of the equation to save.
  Defaults to saving all indices. For example, if you are solving a 3-dimensional ODE, 
  and given `save_idxs = [1, 3]`, only the first and third components of the solution will be outputted. 
  Notice that of course in this case the outputed solution will be two-dimensional.
  
* `tstops`: Denotes *extra* times that the timestepping algorithm must step to.
  This should be used to help the solver deal with discontinuities and
  singularities, since stepping exactly at the time of the discontinuity will
  improve accuracy. If a method cannot change timesteps (fixed timestep multistep methods),
  then `tstops` will use an interpolation, matching the behavior of `saveat`.
  If a method cannot change timesteps and also cannot interpolate,
  then `tstops` must be a multiple of `dt` or else an error will be thrown. Default is `[]`.
* `d_discontinuities:` Denotes locations of discontinuities in low order derivatives.
  This will force FSAL algorithms which assume derivative continuity to re-evaluate
  the derivatives at the point of discontinuity. The default is `[]`.
* `save_everystep`: Saves the result at every timeseries_steps iteration.    
  Default is true if `isempty(saveat)`.
* `timeseries_steps`: Denotes how many steps between saving a value for the
  timeseries. These "steps" are the steps that the solver stops internally (the ones you get by `save_everystep = true`), not the ones that are 
  instructed by the user (all solvers work in a step-like manner). Defaults to 1.
* `save_start`: Denotes whether the initial condition should be included in
  the solution type as the first timepoint. Defaults to true.


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
* `unstable_check`: Specifies a function `unstable_check(dt,t,u)` where, when
  it returns true, it will cause the solver to exit and throw a warning. Defaults
  to `any(isnan,u)`, i.e. checking if any value is a NaN.
* `verbose`: Toggles whether warnings are thrown when the solver exits early.
  Defualts to true.
* `calck`: Turns on and off the internal ability for intermediate    
  interpolations (also known as intermediate density). Not the same as `dense`, which is post-solution interpolation.
  This defaults to `dense || !isempty(saveat) ||  "no custom callback is given"`. 
  This can be used to turn off interpolations
  (to save memory) if one isn't using interpolations when a custom callback is
  used. Another case where this may be used is to turn on interpolations for
  usage in the integrator interface even when interpolations are used nowhere else. 
  Note that this is only required if the algorithm doesn't have
  a free or lazy interpolation (`DP8()`). If `calck = false`, `saveat` cannot be used.
  The rare keyword `calck` can be useful in event handling.

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
