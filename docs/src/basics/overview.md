# Overview of DifferentialEquations.jl

The general workflow for using the package is as follows:

  - Define a problem
  - Solve the problem
  - Analyze the output

## Defining Problems

Problems are specified via a type interface. The problem types are designed to
contain the necessary information to fully define their associated differential
equation. Each problem type has a page explaining their problem type and the special
features associated with them. For example, an ordinary differential equation is defined by

```math
\frac{du}{dt} = f(u,p,t)
```

over some time interval `tspan` with some initial condition `u0`, and therefore
the `ODEProblem` is defined by those components:

```julia
prob = ODEProblem(f, u0, tspan)
prob = ODEProblem(f, u0, tspan, p)
```

Note that the number types in the solution will match the types you designate
in the problem. For example, if one uses `Rational{BigInt}` for specifying the
timespan and `BigFloat` for specifying the initial condition, then the solution
will solve using `Rational{BigInt}` for the timesteps and `BigFloat` for the
independent variables. A wide variety of number types are compatible with the
solvers such as complex numbers, unitful numbers (via Unitful.jl),
decimals (via DecFP), dual numbers, and many more which may not have been tested
yet (thanks to the power of multiple dispatch!). For information on type-compatibility,
please see the solver pages for the specific problems.

## Solving the Problems

Each type of differential equation has its own problem type which allow the solvers
to dispatch to the right methods. The common interface for calling the solvers is:

```julia
sol = solve(prob, alg; kwargs)
```

Into the command, one passes the differential equation problem that they defined
`prob`, optionally choose an algorithm `alg` (a default is given if not
chosen), and change the properties of the solver using keyword arguments. The common
arguments which are accepted by most methods is defined in [the common solver options manual page](@ref solver_options).
The solver returns a solution object `sol` which hold all the details for the solution.

## Analyzing the Solution

With the solution object, you do the analysis as you please! The solution type
has a common interface, which makes handling the solution similar between the
different types of differential equations. Tools such as interpolations
are seamlessly built into the solution interface to make analysis easy. This
interface is described in the [solution handling manual page](@ref solution).

Plotting functionality is provided by a recipe to Plots.jl. To
use plot solutions, simply call the `plot(sol)` and the plotter will generate
appropriate plots. If `save_everystep` was used, the plotters can
generate animations of the solutions to evolution equations using the `animate(sol)`
command. Plots can be customized using all the keyword arguments
provided by Plots.jl. Please see Plots.jl's documentation for more information.

## Add-on Tools

One of the most compelling features of DifferentialEquations.jl is that the
common solver interface allows one to build tools which are “algorithm and
problem agnostic”. For example, one of the provided tools allows for performing
parameter estimation on `ODEProblem`s. Since the `solve` interface is the
same for the different algorithms, one can use any of the associated solving algorithms.
This modular structure allows one to mix and match overarching analysis tools
with specialized algorithms to one's problem, leading to high performance
with a large feature base. Isn't that the promise of Julia just being
fulfilled?

## Development and Testing Tools

Lastly, one unique feature of DifferentialEquations.jl is the existence of algorithm
development and testing functionality. This suite was designed by researchers in
the field of numerical differential equations to both try out new ideas and distribute
finalized results to large audiences. The tools for algorithm development allow for
easy convergence testing, benchmarking, and higher order analysis (stability plotting,
etc.). This is one of the reasons why DifferentialEquations.jl contains many algorithms
which are unique and the results of recent publications! Please check out the
[developer documentation](https://devdocs.sciml.ai/dev/)
for more information on using the development tools.

Note that DifferentialEquations.jl allows for distributed development, meaning that
algorithms which “plug-into the ecosystem” don't have to be a part of the major packages.
If you are interested in adding your work to the ecosystem, checkout the [developer documentation](https://devdocs.sciml.ai/dev/)
for more information.
