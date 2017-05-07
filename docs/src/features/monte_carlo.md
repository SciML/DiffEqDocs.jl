# Parallel Monte Carlo Simulations

## Performing a Monte Carlo Simulation

### Building a Problem

To perform a Monte Carlo simulation, define a `MonteCarloProblem`. The constructor is:

```julia
MonteCarloProblem(prob::DEProblem;
                  output_func = (sol,i) -> sol,
                  prob_func= (prob,i)->prob)
```

* `prob_func`: The function by which the problem is to be modified.
* `output_func`: The reduction function.

One can specify a function `prob_func` which changes the problem. For example:

```julia
function prob_func(prob,i)
  prob.u0 = randn()*prob.u0
  prob
end
```

modifies the initial condition for all of the problems by a standard normal
random number (a different random number per simulation). This can be used
to perform searches over initial values. Note that the parameter `i` is a unique
counter over the simulations. Thus if you have an array of initial conditions `u0_arr`,
you can have the `i`th simulation use the `i`th initial condition via:

```julia
function prob_func(prob,i)
  prob.u0 = u0_arr[i]
  prob
end
```

If your function is a `ParameterizedFunction`,
you can do similar modifications to `prob.f` to perform a parameter search. The `output_func`
is a reduction function. It's arguments are the generated solution and the unique
index for the run. For example, if we wish to only save the 2nd coordinate
at the end of each solution, we can do:

```julia
output_func(sol,i) = sol[end,2]
```

Thus the Monte Carlo Simulation would return as its data an array which is the
end value of the 2nd dependent variable for each of the runs.

### Parameterizing the Monte Carlo Components

The Monte Carlo components can be parameterized by using the `ConcreteParameterizedFunction`
constructors.

```julia
ProbParameterizedFunction(prob_func,params)
OutputParameterizedFunction(output_func,params)
```

Here, the signatures are `prob_func(prob,i,params)` and `output_func(sol,params)`.
These parameters are added to the parameter list for use in the parameter estimation
schemes.

### Solving the Problem

```julia
sim = solve(prob,alg,kwargs...)
```

The keyword arguments take in the arguments for the common solver interface and will
pass them to the differential equation solver. The special keyword arguments to note are:

* `num_monte`: The number of simulations to run. Default is 10,000.
* `parallel_type` : The type of parallelism to employ. Default is `:pmap`.

The types of parallelism included are:

* `:none` - No parallelism
* `:threads` - This uses multithreading. It's local (single computer, shared memory)
  parallelism only. Fastest when the trajectories are quick.
* `:parfor` - A multiprocessing parallelism. Slightly better than `pmap` when the
  calculations are fast. Does not re-distribute work: each trajectory is assumed
  to take as long to calculate.
* `:pmap` - The default. Uses `pmap` internally. It will use as many processors as you
  have Julia processes. To add more processes, use `addprocs(n)`. See Julia's
  documentation for more details. Recommended for the case when each trajectory
  calculation isn't "too quick" (at least about a millisecond each?).
* `:split_threads` - This uses threading on each process, splitting the problem
  into `nprocs()` even parts. This is for solving many quick trajectories on a
  multi-node machine. It's recommended you have one process on each node.

Additionally, a `MonteCarloEstimator` can be supplied

```julia
sim = solve(prob,estimator,alg,kwargs...)
```

These will be detailed when implemented.

### Solution Type

The resulting type is a `MonteCarloSimulation`, which includes the array of
solutions. If the problem was a `TestProblem`, summary statistics on the errors
are returned as well.

### Plot Recipe

There is a plot recipe for a `AbstractMonteCarloSimulation` which composes all
of the plot recipes for the component solutions. The keyword arguments are passed
along. A useful argument to use is `linealpha` which will change the transparency
of the plots. An additional argument is `idxs` which allows you to choose which
components of the solution to plot. For example, if the differential equation
is a vector of 9 values, `idxs=1:2:9` will plot only the Monte Carlo solutions
of the odd components.

### Example

Let's test the sensitivity of the linear ODE to its initial condition.

```julia
addprocs(4)
using DiffEqMonteCarlo, DiffEqBase, DiffEqProblemLibrary, OrdinaryDiffEq
prob = prob_ode_linear
prob_func = function (prob)
  prob.u0 = rand()*prob.u0
  prob
end
monte_prob = MonteCarloProblem(prob,prob_func=prob_func)
sim = solve(monte_prob,Tsit5(),num_monte=100)

using Plots
plotly()
plot(sim,linealpha=0.4)
```

Here we solve the same ODE 100 times on 4 different cores, jiggling the initial
condition by `rand()`. The resulting plot is as follows:

![monte_carlo_plot](../assets/monte_carlo_plot.png)

## Analyzing a Monte Carlo Experiment

Analysis tools are included for generating summary statistics and summary plots
for a `MonteCarloSimulation`.

### Time steps vs time points

For the summary statistics, there are two types. You can either summarize by
time steps or by time points. Summarizing by time steps assumes that the time steps
are all the same time point, i.e. the integrator used a fixed `dt` or the values were
saved using `saveat`. Summarizing by time points requires interpolating the solution.

### Analysis at a time step or time point

```julia
get_timestep(sim,i) # Returns an iterator of each simulation at time step i
get_timepoint(sim,t) # Returns an iterator of each simulation at time point t
componentwise_vectors_timestep(sim,i) # Returns a vector of each simulation at time step i
componentwise_vectors_timepoint(sim,t) # Returns a vector of each simulation at time point t
```

### Summary Statistics

#### Single Time Statistics

The available functions for time steps are:

```julia
timestep_mean(sim,i) # Computes the mean of each component at time step i
timestep_median(sim,i) # Computes the median of each component at time step i
timestep_quantile(sim,q,i) # Computes the quantile q of each component at time step i
timestep_meanvar(sim,i)  # Computes the mean and variance of each component at time step i
timestep_meancov(sim,i,j) # Computes the mean at i and j, and the covariance, for each component
timestep_meancor(sim,i,j) # Computes the mean at i and j, and the correlation, for each component
timestep_weighted_meancov(sim,W,i,j) # Computes the mean at i and j, and the weighted covariance W, for each component
```

The available functions for time points are:

```julia
timepoint_mean(sim,t) # Computes the mean of each component at time t
timepoint_median(sim,t) # Computes the median of each component at time t
timepoint_quantile(sim,q,t) # Computes the quantile q of each component at time t
timepoint_meanvar(sim,t) # Computes the mean and variance of each component at time t
timepoint_meancov(sim,t1,t2) # Computes the mean at t1 and t2, the covariance, for each component
timepoint_meancor(sim,t1,t2) # Computes the mean at t1 and t2, the correlation, for each component
timepoint_weighted_meancov(sim,W,t1,t2) # Computes the mean at t1 and t2, the weighted covariance W, for each component
```

#### Full Timeseries Statistics

Additionally, the following functions are provided for analyzing the full timeseries.
The `mean` and `meanvar` versions return a `DiffEqArray` which can be directly plotted.
The `meancov` and `meancor` return a matrix of tuples, where the tuples are the
`(mean_t1,mean_t2,cov or cor)`.

The available functions for the time steps are:

```julia
timeseries_steps_mean(sim) # Computes the mean at each time step
timeseries_steps_median(sim) # Computes the median at each time step
timeseries_steps_quantile(sim,q) # Computes the quantile q at each time step
timeseries_steps_meanvar(sim) # Computes the mean and variance at each time step
timeseries_steps_meancov(sim) # Computes the covariance matrix and means at each time step
timeseries_steps_meancor(sim) # Computes the correlation matrix and means at each time step
timeseries_steps_weighted_meancov(sim) # Computes the weighted covariance matrix and means at each time step
```

The available functions for the time points are:

```julia
timeseries_point_mean(sim,ts) # Computes the mean at each time point in ts
timeseries_point_median(sim,ts) # Computes the median at each time point in ts
timeseries_point_quantile(sim,q,ts) # Computes the quantile q at each time point in ts
timeseries_point_meanvar(sim,ts) # Computes the mean and variance at each time point in ts
timeseries_point_meancov(sim,ts) # Computes the covariance matrix and means at each time point in ts
timeseries_point_meancor(sim,ts) # Computes the correlation matrix and means at each time point in ts
timeseries_point_weighted_meancov(sim,ts) # Computes the weighted covariance matrix and means at each time point in ts
```

### MonteCarloSummary

The `MonteCarloSummary` type is included to help with analyzing the general summary
statistics. Two constructors are provided:

```julia
MonteCarloSummary(sim)
MonteCarloSummary(sim,ts)
```

The first produces a `(mean,var)` summary at each time step. As with the summary
statistics, this assumes that the time steps are all the same. The second produces
a `(mean,var)` summary at each time point `t` in `ts`. This requires the ability
to interpolate the solution.

#### Plot Recipe

The `MonteCarloSummary` comes with a plot recipe for visualizing the summary
statistics. The extra keyword arguments are:

- `idxs`: the solution components to plot. Defaults to plotting all components.
- `error_style`: The style for plotting the error. Defaults to `ribbon`. Other
  choices are `:bars` for error bars and `:none` for no error bars.

One useful argument is `fillalpha` which controls the transparency of the ribbon
around the mean. The confidence interval is the Gaussian CI `1.96*var`.

### Example Analysis

In this example we will show how to analyze a `MonteCarloSolution`. First, let's
generate a 10 solution Monte Carlo experiment:

```julia
prob = prob_sde_2Dlinear
prob2 = MonteCarloProblem(prob)
sim = solve(prob2,SRIW1(),dt=1//2^(3),num_monte=10,adaptive=false)
```

The system, `prob_sde_2Dlinear`, is a `(4,2)` system of stochastic
differential equations which we solved 10 times. We can compute the
mean and the variance at the 3rd timestep using:

```julia
m,v = timestep_meanvar(sim,3)
```

or we can compute the mean and the variance at the `t=0.5` using:

```julia
m,v = timepoint_meanvar(sim,0.5)
```

We can get a series for the mean and the variance at each time step using:

```julia
m_series,v_series = timeseries_steps_meanvar(sim)
```

or at chosen values of `t`:

```julia
ts = 0:0.1:1
m_series = timeseries_point_mean(sim,ts)
```

Note that these mean and variance series can be directly plotted. We can
compute covariance matrices similarly:

```julia
timeseries_steps_meancov(sim) # Use the time steps, assume fixed dt
timeseries_point_meancov(sim,0:1//2^(3):1,0:1//2^(3):1) # Use time points, interpolate
```

For general analysis, we can build a `MonteCarloSummary` type.

```julia
summ = MonteCarloSummary(sim)
```

will summarize at each time step, while

```julia
summ = MonteCarloSummary(sim,0.0:0.1:1.0)
```

will summarize at the `0.1` time points using the interpolations. To
visualize the results we can plot it. Since there are 8 components to
the differential equation, this can get messy, so let's only plot the
3rd component:

```julia
plot(summ;idxs=3)
```

![monte_ribbon](../assets/monte_ribbon.png)

We can change to errorbars instead of ribbons and plot two different
indices:

```julia
plot(summ;idxs=(3,5),error_style=:bars)
```

![monte_bars](../assets/monte_bars.png)

Or we can simply plot the mean of every component over time:

```julia
plot(summ;error_style=:none)
```

![monte_means](../assets/monte_means.png)
