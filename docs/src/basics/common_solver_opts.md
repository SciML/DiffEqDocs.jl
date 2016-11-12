# Common Solver Options

The DifferentialEquations.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. For more detailed information on the defaults and the available options
for specific algorithms / packages, see the [Option Availability manual page]().

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
way to only saving the solution at the final timepoint. For more information on
controlling the output behavior, see the [Output Specification manual page](../man/output_specification.html).

* `dense`: Denotes whether to save the extra pieces for dense (continuous) output. Default is true
  for algorithms which have the ability to produce dense output.
* `saveat`: Denotes extra times to save the solution at during the solving phase. Note that this
  can be used even if `dense=false`. Default is `[]`.
* `tstops`: Denotes extra times that the timestepping algorithm must step to. This should
  only be used if dense output via `saveat` is not available for the algorithm (for efficiency).
  Default is `[]`.
* `calck`: Turns on and off the internal ability for intermediate interpolations. This defaults
  to `dense || !isempty(saveat) || `"no custom callback is given". This can be used
  to turn off interpolations (to save memory) even when a custom callback is used.
* `save_timeseries`: Saves the result at every timeseries_steps iteration. Default is true.
* `timeseries_steps`: Denotes how many steps between saving a value for the timeseries. Defaults to 1.

## Stepsize Control

These arguments control the timestepping routines.

* `adaptive` - Turns on adaptive timestepping for appropriate methods. Default is true.
* `abstol` - Absolute tolerance in adaptive timestepping. Defaults to 1e-3.
* `reltol` - Relative tolerance in adaptive timestepping. Defaults to 1e-6.
* `dt`: Sets the initial stepsize. This is also the stepsize for fixed timestep methods.
  Defaults to an automatic choice.
* `internalnorm` - The norm function `internalnorm(u)` which error estimates are calculated.
  Defaults are package-dependent.
* `gamma` - The risk-factor γ in the q equation for adaptive timestepping. Default is algorithm dependent.
* `dtmax` - Maximum dt for adaptive timestepping. Defaults are package-dependent.
* `dtmin` - Minimum dt for adaptive timestepping. Defaults are package-dependent.
* `beta1` - The Lund stabilization α parameter. Defaults are algorithm-dependent.
* `beta2` - The Lund stabilization β parameter. Defaults are algorithm-dependent.
* `qmax` - Defines the maximum value possible for the adaptive q. Defaults are algorithm-dependent.
* `qmin` - Defines the maximum value possible for the adaptive q. Defaults are algorithm-dependent.
* `qoldinit` - The initial `qold` in stabilization stepping. Defaults are algorithm-dependent.

## Miscellaneous

* `maxiters` - Maximum number of iterations before stopping. Defaults to 1e5.
* `autodiff` - Turns on/off the use of autodifferentiation (via ForwardDiff) in the
  implicit solvers which use `NLsolve`. Default is true.
* `callback` - Specifies a callback function. Defaults to a callback function which
  performs the saving routine. For more information, see the
  [Event Handling and Callback Functions manual page](https://juliadiffeq.github.io/DiffEqDocs.jl/latest/man/callback_functions.html).

## Progress Bar Control

These arguments control the usage of the progressbar in the Juno IDE.

* `progressbar` - Turns on/off the Juno progressbar. Default is false.
* `progress_steps` - Numbers of steps between updates of the progress bar. Default is 1000.
* `progressbar_name` - Controls the name of the progressbar. Default is the name of the problem type.

## Error Calculations

If you are using the test problems (ex: `ODETestProblem`), then the following options
control the errors which are calculated:

* `timeseries_errors` - Turns on and off the calculation of errors at the steps which
  were taken, such as the `l2` error. Default is true.
* `dense_errors` - Turns on and off the calculation of errors at the steps which
  require dense output and calculate the error at 100 evenly-spaced points throughout
  `tspan`. An example is the `L2` error. Default is false.
