# DifferentialEquations.jl Documentation

DifferentialEquations.jl is a package for numerically solving differential
equations using the various components of JuliaDiffEq. The purpose of this package is to supply
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

If you have any questions, or just want to chat about solvers/using the package,
please feel free to use the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby).
For bug reports, feature requests, etc., please submit an issue. If you're
interested in contributing, please see the
[Developer Documentation](https://juliadiffeq.github.io/DiffEqDevDocs.jl/latest/).

### Using the Package

To install the package, use the following command inside the Julia REPL:
```julia
Pkg.add("DifferentialEquations")
```

To load the package, use the command:

```julia
using DifferentialEquations
```

To understand the package in more detail, check out the following tutorials in the manual.
Example IJulia notebooks [can be found in DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl).
If you find any example where there seems to be an error, please open an issue.

For the most up to date on using the package information, please join [the Gitter channel](https://gitter.im/JuliaDiffEq/Lobby)

Using the bleeding edge for the latest features and development is only recommended
for power users. Information for how to get to the bleeding edge is found in the
[developer documentation](https://juliadiffeq.github.io/DiffEqDevDocs.jl/latest/index.html#Bleeding-Edge-1).

### Supported Equations

For PDEs, one can optionally specify a noise equation. The solvers currently have
stochastic variants for handling Gaussian Space-time white noise SPDEs.

* Ordinary Differential Equations (ODEs)
* Stochastic Differential Equations (SDEs)
* Algebraic Differential Equations (DAEs)
* Delay Differential Equations (DDEs)
* (Stochastic) Partial Differential Equations ((S)PDEs):

  * Linear Poisson Equation
  * Semi-linear Poisson Equation
  * Linear Heat Equation
  * Semi-linear Heat Equation (aka Reaction-Diffusion Equation)
  * Stationary Stokes Equation

For help with choosing a solver algorithm, please see the solver options pages.

### IJulia Notebook Tutorials

You can access extra tutorials supplied in the [DiffEqTutorials.jl repository](https://github.com/JuliaDiffEq/DiffEqTutorials.jl).
If you have [IJulia](https://github.com/JuliaLang/IJulia.jl) installed, you can
view them locally and interactively, by cloning the repository:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqTutorials.jl")
using IJulia
notebook(dir = Pkg.dir("DiffEqTutorials"))
```

### Tutorials

The following tutorials will introduce you to the functionality of DifferentialEquations.jl
More examples can be found by [checking out the IJulia notebooks in the examples folder](https://github.com/JuliaDiffEq/DiffEqTutorials.jl).

```@contents
Pages = [
    "tutorials/ode_example.md",
    "tutorials/sde_example.md",
    "tutorials/dde_example.md",
    "tutorials/dae_example.md",
    "tutorials/fempoisson_example.md",
    "tutorials/femheat_example.md",
    "tutorials/femstochastic_example.md"
    ]
Depth = 2
```

### Basics

These pages introduce you to the core of DifferentialEquations.jl and the common
interface. It explains the general workflow, options which are generally available,
and the general tools for analysis.

```@contents
Pages = [
    "basics/overview.md",
    "basics/common_solver_opts.md",
    "basics/solution.md",
    "basics/plot.md",
    "basics/integrator.md"
    ]
Depth = 2
```

### Equation Types

These pages describe building the problem types to define differential equations
for the solvers, and the special features of the different solution types.

```@contents
Pages = [
  "types/ode_types.md",
  "types/sde_types.md",
  "types/dde_types.md",
  "types/dae_types.md",
  "types/fem_types.md",
  "types/stokes_types.md"
]
Depth = 2
```

### Solver Algorithms

These pages describe the solvers and available algorithms in detail.

```@contents
Pages = [
  "solvers/ode_solve.md",
  "solvers/sde_solve.md",
  "solvers/dde_solve.md",
  "solvers/dae_solve.md",
  "solvers/fempoisson_solve.md",
  "solvers/femheat_solve.md",
  "solvers/fdmstokes_solve.md"
]
Depth = 2
```

### Additional Features

These sections discuss extra performance enhancements, event handling, and other
in-depth features.

```@contents
Pages = [
    "features/performance_overloads.md",
    "features/callback_functions.md",
    "features/callback_library.md",
    "features/monte_carlo.md",
    "features/mesh.md",
    "features/output_specification.md",
    "features/progress_bar.md"
]
Depth = 2
```

### Analysis Tools

Because DifferentialEquations.jl has a common interface on the solutions, it is
easy to add functionality to the entire DiffEq ecosystem by developing it
to the solution interface. These pages describe the add-on analysis tools which
are available.

```@contents
Pages = [
    "analysis/parameterized_functions.md",
    "analysis/parameter_estimation.md",
    "analysis/sensitivity.md",
    "analysis/uncertainty_quantification.md",
    "analysis/dev_and_test.md"
]
Depth = 2
```
