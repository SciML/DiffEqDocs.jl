# DifferentialEquations.jl Documentation

DifferentialEquations.jl is a package for solving numerically solving differential
equations. The purpose of this package is to supply
efficient Julia implementations of solvers for various differential equations.
Equations within the realm of this package include ordinary differential equations
(ODEs), stochastic ordinary differential equations (SODEs or SDEs), stochastic
partial differential equations (SPDEs), partial differential equations (with both
finite difference and finite element methods), differential algebraic equations
(DAEs), and differential delay equations (DDEs). The well-optimized
DifferentialEquations solvers benchmark as the fastest Julia implementations,
using classic algorithms and ones from recent research, and include algorithms
optimized for high-precision and HPC applications.  It integrates with the
Julia package sphere, for example using Juno's progress meter, automatic plotting,
built-in interpolations, and wraps other differential equation solvers so
that many different methods for solving the equations can be accessed by
simply switching a keyword argument. It utilizes Julia's generality to be
able to solve problems specified with arbitrary number types (types with
units like Unitful, and arbitrary precision numbers like BigFloats and
ArbFloats), arbitrary sized arrays (ODEs on matrices), and more. This gives
a powerful mixture of speed and productivity features to help you solve and
analyze your differential equations faster.

Example IJulia notebooks [can be found in the examples folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).
If you find any example where there seems to be an error, please open an issue.

If you have any questions, or just want to chat about solvers/using the package, please feel free to use the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby). For bug reports, feature requests, etc., please submit an issue. If you're interested in contributing, please see the [Contributor's Guide](http://juliadiffeq.github.io/DifferentialEquations.jl/latest/internals/contributors_guide.html).

## Using the Package

To install the package, use the following command inside the Julia REPL:
```julia
Pkg.add("DifferentialEquations")
```

To load the package, use the command:

```julia
using DifferentialEquations
```

To understand the package in more detail, check out the following tutorials in the manual. Examples
IJulia notebooks using DifferentialEquations can be found [in the examples folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).
Codes for the latest features can be found in [test/](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/test).

For the most up to date on using the package information, please contact me [via the repository Gitter](https://gitter.im/JuliaDiffEq/Lobby)
or [read the latest documentation](http://JuliaDiffEq.github.io/DifferentialEquations.jl/latest/)

## Bleeding Edge

This package suite is still under heavy development. If you are a power user
and would like to try out the latest features, it is recommended you use the
MetaDiffEq metapackage. To do so, use the following commands:

```julia
Pkg.clone("https://github.com/tbreloff/MetaPkg.jl") # Install MetaPkg
using MetaPkg
meta_add("MetaDiffEq") # Adds all of the packages, even those unregistered
meta_checkout("MetaDiffEq") # Checks out the master branch on all of the packages
```

Note that this is for power users who are familiar with Julia. If you are having
issues, please contact Chris Rackauckas in  [the Gitter channel.](https://gitter.im/JuliaDiffEq/Lobby)

## Supported Equations

For PDEs, one can optionally specify a noise equation. The solvers currently have
stochastic variants for handling Gaussian Space-time white noise SPDEs.

* ODEs
* SODEs
* DAEs
* (Stochastic) PDEs

  * Linear Poisson Equation
  * Semi-linear Poisson Equation
  * Linear Heat Equation
  * Semi-linear Heat Equation (aka Reaction-Diffusion Equation)
  * Stationary Stokes Equation

For help with choosing a solver algorithm, please see the solver options pages.

## IJulia Notebook Tutorials

You can access extra tutorials supplied in the [DiffEqTutorials.jl repository](https://github.com/JuliaDiffEq/DiffEqTutorials.jl).
If you have [IJulia](https://github.com/JuliaLang/IJulia.jl) installed, you can
view them locally and interactively, by cloning the repository:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqTutorials.jl")
using IJulia
notebook(dir = Pkg.dir("DiffEqTutorials"))
```

## Tutorials

The following tutorials will introduce you to the functionality of DifferentialEquations.jl
More examples can be found by [checking out the IJulia notebooks in the examples
folder](https://github.com/JuliaDiffEq/DifferentialEquations.jl/tree/master/examples).

```@contents
Pages = [
    "tutorials/ode_example.md",
    "tutorials/sde_example.md",
    "tutorials/dae_example.md",
    "tutorials/fempoisson_example.md",
    "tutorials/femheat_example.md",
    "tutorials/femstochastic_example.md"
    ]
Depth = 2
```

## Solver Options

These pages describe the options available in the solvers.

```@contents
Pages = [
  "solvers/common_solvers_opts.md"
  "solvers/ode_solve.md",
  "solvers/sde_solve.md",
  "solvers/dae_solve.md",
  "solvers/fempoisson_solve.md",
  "solvers/femheat_solve.md",
  "solvers/fdmstokes_solve.md"
]
Depth = 2
```

## Manual

```@contents
Pages = [
    "man/overview.md",
    "man/ODEProblem.md",
    "man/SDEProblem.md",
    "man/FEMProblem.md",
    "man/StokesProblem.md",
    "man/mesh.md",
    "man/solution.md",
    "man/output_specification.md",
    "man/callback_functions.md",
    "man/plot.md",
    "man/parameter_estimation.md",
    "man/sensitivity.md",
    "man/function_definition_macros.md",
    "man/benchmarks.md",
    "man/convergence.md",
    "man/conditional_dependencies.md",
    "man/progress_bar.md"
]
Depth = 2
```

## Internal Documentation

```@contents
Pages = [
  "internals/contributors_guide.md",
  "internals/fem_tools.md",
  "internals/extras.md",
  "internals/solver_helpers.md",
  "internals/notes_on_algorithms.md",
  "internals/function_index.md"
]
Depth = 2
```
