# Parallel Monte Carlo Simulations

DiffEqMonteCarlo.jl provides functionality for easily performing parallel Monte
Carlo simulations using the DiffEq solvers. As the package is currently unregistered,
it must first be installed via:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqMonteCarlo.jl")
using DiffEqMonteCarlo
```

## Performing a Monte Carlo Simulation

To perform a Monte Carlo simulation, you simply use the interface:

```julia
sim = monte_carlo_simulation(prob,alg,kwargs...)
```

The keyword arguments take in the arguments for the common solver interface.
The special keyword arguments to note are:

* num_monte: The number of simulations to run

In addition, one can specify a function `u0_func` which changes the initial
condition around. For example:

```julia
function prob_func(prob)
  prob.u0 = randn()*prob.u0
  prob
end
```

Modifies the initial condition for all of the problems by a standard normal
random number (a different random number per simulation). This can be used
to perform searches over initial values. If your function is a `ParameterizedFunction`,
you can do similar modifications to `f` to perform a parameter search. One then
passes this function via:

```julia
sim = monte_carlo_simulation(prob,alg,prob_func,kwargs...)
```

## Parallelism

Since this is using `pmap` internally, it will use as many processors as you
have Julia processes. To add more processes, use `addprocs(n)`. See Julia's
documentation for more details.

## Solution

The resulting type is a `MonteCarloSimulation`, which includes the array of
solutions. If the problem was a `TestProblem`, summary statistics on the errors
are returned as well.

## Plot Recipe

There is a plot recipe for a `AbstractMonteCarloSimulation` which composes all
of the plot recipes for the component solutions. The keyword arguments are passed
along. A useful argument to use is `linealpha` which will change the transparency
of the plots.

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
sim = monte_carlo_simulation(prob,Tsit5(),prob_func,num_monte=100)

using Plots
plotly()
plot(sim,linealpha=0.4)
```

Here we solve the same ODE 100 times on 4 different cores, jiggling the initial
condition by `rand()`. The resulting plot is as follows:

![monte_carlo_plot](../assets/monte_carlo_plot.png)
