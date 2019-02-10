# Common Solver Options

The DifferentialEquations.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. Not all of the interface is provided by every algorithm.
For more detailed information on the defaults and the available options
for specific algorithms / packages, see the manual pages for the solvers of specific
problems. To see whether a specific package is compaible with the use of a
given option, see the [Solver Compatibility Chart](@ref)

## Default Algorithm Hinting

To help choose the default algorithm, the keyword argument `alg_hints` is
provided to `solve`. `alg_hints` is a `Vector{Symbol}` which describe the
problem at a high level to the solver. The options are:

* `:auto` vs `:nonstiff` vs `:stiff` - Denotes the equation as nonstiff/stiff.
  `:auto` allow the default handling algorithm to choose stiffness detection
  algorithms. The default handling defaults to using `:auto`.

Currently unused options include:

* `:interpolant` - Denotes that a high-precision interpolation is important.
* `:memorybound` - Denotes that the solver will be memory bound.

This functionality is derived via the benchmarks in
[DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl)

### SDE Specific Alghints

* `:additive` - Denotes that the underlying SDE has additive noise.
* `:stratonovich` - Denotes that the solution should adhere to the Stratonovich
  interpretation.

## Output Control

These arguments control the output behavior of the solvers. It defaults to maximum
output to give the best interactive user experience, but can be reduced all the
way to only saving the solution at the final timepoint.

The following options are all related to output control. See the "Examples"
section at the end of this page for some example usage.

* `dense`: Denotes whether to save the extra pieces required for dense (continuous)
  output. Default is true for algorithms which have the ability to produce dense
  output. If dense is false, the solution still acts like a function, and `sol(t)`
  is a linear interpolation between the saved time points.
* `saveat`: Denotes specific times to save the solution at, during the solving
  phase. The solver will save at each of the timepoints in this array in the
  most efficient manner (always including the points of `tspan`). Note that this
  can be used even if `dense=false`. In fact, if only `saveat` is given, then
  the arguments `save_everystep` and `dense` are becoming `false` by default and
  must be explicitly given as `true` if desired.   If `saveat` is given a number,
  then it will automatically expand to `tspan[1]:saveat:tspan[2]`. For methods
  where interpolation is not possible, `saveat` may be equivalent to `tstops`.
  Default is `[]`.
* `save_idxs`: Denotes the indices for the components of the equation to save.
  Defaults to saving all indices. For example, if you are solving a 3-dimensional ODE,
  and given `save_idxs = [1, 3]`, only the first and third components of the
  solution will be outputted.
  Notice that of course in this case the outputed solution will be two-dimensional.
* `tstops`: Denotes *extra* times that the timestepping algorithm must step to.
  This should be used to help the solver deal with discontinuities and
  singularities, since stepping exactly at the time of the discontinuity will
  improve accuracy. If a method cannot change timesteps (fixed timestep
  multistep methods), then `tstops` will use an interpolation,
  matching the behavior of `saveat`. If a method cannot change timesteps and
  also cannot interpolate, then `tstops` must be a multiple of `dt` or else an
  error will be thrown. Default is `[]`.
* `d_discontinuities:` Denotes locations of discontinuities in low order derivatives.
  This will force FSAL algorithms which assume derivative continuity to re-evaluate
  the derivatives at the point of discontinuity. The default is `[]`.
* `save_everystep`: Saves the result at every step.    
  Default is true if `isempty(saveat)`.
* `save_on`: Denotes whether intermediate solutions are saved. This overrides the
  settings of `dense`, `saveat` and `save_everystep` and is used by some applicatioins
  to manually turn off saving temporarily. Everyday use of the solvers should leave
  this unchanged. Defaults to `true`.
* `save_start`: Denotes whether the initial condition should be included in
  the solution type as the first timepoint. Defaults to `true`.
* `save_end`: Denotes whether the final timepoint is forced to be saved,
  regardless of the other saving settings. Defaults to `true`.
* `initialize_save`: Denotes whether to save after the callback initialization
  phase (when `u_modified=true`). Defaults to `true`.

Note that `dense` requires `save_everystep=true` and `saveat=false`. If you need
additional saving while keeping dense output, see
[the SavingCallback in the Callback Library](http://docs.juliadiffeq.org/latest/features/callback_library.html).

## Stepsize Control

These arguments control the timestepping routines.

#### Basic Stepsize Control

These are the standard options for controlling stepping behavior. Error estimates
do the comparison

```math
err_{scaled} = err/(abstol + max(uprev,u)*reltol)
```

The scaled error is guaranteed to be `<1` for a given local error estimate
(note: error estimates are local unless the method specifies otherwise). `abstol`
controls the non-scaling error and thus can be though of as the error around zero.
`reltol` scales with the size of the dependent variables and so one can interpret
`reltol=1e-3` as roughly being (locally) correct to 3 digits. Note tolerances can
be specified element-wise by passing a vector whose size matches `u0`.

* `adaptive`: Turns on adaptive timestepping for appropriate methods. Default
  is true.
* `abstol`: Absolute tolerance in adaptive timestepping. This is the tolerance
  on local error estimatoes, not necessarily the global error (though these quantities
  are related). Defaults to `1e-6` on deterministic equations (ODEs/DDEs/DAEs) and `1e-2`
  on stochastic equations (SDEs/RODEs).
* `reltol`: Relative tolerance in adaptive timestepping.  This is the tolerance
  on local error estimatoes, not necessarily the global error (though these quantities
  are related). Defaults to `1e-3` on deterministic equations (ODEs/DDEs/DAEs) and `1e-2`
  on stochastic equations (SDEs/RODEs).
* `dt`: Sets the initial stepsize. This is also the stepsize for fixed
  timestep methods. Defaults to an automatic choice if the method is adaptive.
* `dtmax`: Maximum dt for adaptive timestepping. Defaults are
  package-dependent.
* `dtmin`: Minimum dt for adaptive timestepping. Defaults are
  package-dependent.
* `force_dtmin`: Declares whether to continue, forcing the minimum `dt` usage.
  Default is `false`, which has the solver throw a warning and exit early when
  encountering the minimum `dt`. Setting this true allows the solver to continue,
  never letting `dt` go below `dtmin` (and ignoring error tolerances in those
  cases). Note that `true` is not compatible with most interop packages.

#### Fixed Stepsize Usage

Note that if a method does not have adaptivity, the following rules apply:

* If `dt` is set, then the algorithm will step with size `dt` each iteration.
* If `tstops` and `dt` are both set, then the algorithm will step with either a
  size `dt`, or use a smaller step to hit the `tstops` point.
* If `tstops` is set without `dt`, then the algorithm will step directly to
  each value in `tstops`
* If neither `dt` nor `tstops` are set, the solver will throw an error.

#### Advanced Adaptive Stepsize Control

These arguments control more advanced parts of the internals of adaptive timestepping
and are mostly used to make it more efficient on specific problems. For detained
explanations of the timestepping algorithms, see the
[timestepping descriptions](../../extras/timestepping.html)

* `internalnorm`: The norm function `internalnorm(u,t)` which error estimates
  are calculated. Required are two dispatches: one dispatch for the state variable
  and the other on the elements of the state variable (scalar norm).
  Defaults are package-dependent.
* `gamma`: The risk-factor γ in the q equation for adaptive timestepping.
  Default is algorithm dependent.
* `beta1`: The Lund stabilization α parameter. Defaults are
  algorithm-dependent.
* `beta2`: The Lund stabilization β parameter. Defaults are
  algorithm-dependent.
* `qmax`: Defines the maximum value possible for the adaptive q. Defaults are
  algorithm-dependent.
* `qmin`: Defines the maximum value possible for the adaptive q. Defaults are
  algorithm-dependent.
* `qsteady_min`: Defines the minimum for the range around 1 where the timestep
  is held constant. Defaults are algorithm-dependent.
* `qsteady_max`: Defines the maximum for the range around 1 where the timestep
  is held constant. Defaults are algorithm-dependent.
* `qoldinit`: The initial `qold` in stabilization stepping. Defaults are
  algorithm-dependent.
* `failfactor`: The amount to decrease the timestep by if the Newton iterations
  of an implicit method fail. Default is 2.

## Memory Optimizations

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
* `alias_u0`: allows the solver to alias the initial condition array that is contained
  in the problem struct. Defaults to false.

## Miscellaneous

* `maxiters`: Maximum number of iterations before stopping. Defaults to 1e5.
* `callback`: Specifies a callback. Defaults to a callback function which
  performs the saving routine. For more information, see the
  [Event Handling and Callback Functions manual page](../../features/callback_functions.html).
* `isoutofdomain`: Specifies a function `isoutofdomain(u,p,t)` where, when it
  returns true, it will reject the timestep. Disabled by default.
* `unstable_check`: Specifies a function `unstable_check(dt,u,p,t)` where, when
  it returns true, it will cause the solver to exit and throw a warning. Defaults
  to `any(isnan,u)`, i.e. checking if any value is a NaN.
* `verbose`: Toggles whether warnings are thrown when the solver exits early.
  Defualts to true.

## Progress Monitoring

These arguments control the usage of the progressbar in the Juno IDE.

* `progress`: Turns on/off the Juno progressbar. Default is false.
* `progress_steps`: Numbers of steps between updates of the progress bar.
  Default is 1000.
* `progress_name`: Controls the name of the progressbar. Default is the name
  of the problem type.
* `progress_message`: Controls the message with the progressbar. Defaults to
  showing `dt`, `t`, the maximum of `u`.

## Error Calculations

If you are using the test problems (ex: `ODETestProblem`), then the following
options control the errors which are calculated:

* `timeseries_errors`: Turns on and off the calculation of errors at the steps
  which were taken, such as the `l2` error. Default is true.
* `dense_errors`: Turns on and off the calculation of errors at the steps which
  require dense output and calculate the error at 100 evenly-spaced points
  throughout `tspan`. An example is the `L2` error. Default is false.

## Examples

The following lines are examples of how one could use the configuration of
`solve()`. For these examples a 3-dimensional ODE problem is assumed, however
the extention to other types is straightforward.

1. `solve(prob, AlgorithmName())` : The "default" setting, with a user-specified
  algorithm (given by `AlgorithmName()`). All parameters get their default values.
  This means that the solution is saved at the steps the Algorithm stops internally
  and dense output is enabled if the chosen algorithm allows for it.

  All other integration parameters (e.g. stepsize) are chosen automatically.
2. `solve(prob, saveat = 0.01, abstol = 1e-9, reltol = 1e-9)` : Standard setting
  for accurate output at specified (and equidistant) time intervals, used for
  e.g. Fourier Transform. The solution is given every 0.01 time units,
  starting from `tspan[1]`. The solver used is `Tsit5()` since no keyword
  `alg_hits` is given.

3. `solve(prob, maxiters = 1e7, progress = true, save_idxs = [1])` : Using longer
  maximum number of solver iterations can be useful when a given `tspan` is very
  long. This example only saves the first of the variables of the system, either
  to save size or because the user does not care about the others. Finally, with
  `progress = true` you are enabling the progress bar, provided you are using
  the Atom+Juno IDE set-up for your Julia.
