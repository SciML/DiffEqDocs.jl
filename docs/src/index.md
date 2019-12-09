# DifferentialEquations.jl Documentation

This is a suite for numerically solving differential equations written in Julia
and available for use in Julia, Python, and R. The
purpose of this package is to supply efficient Julia implementations of solvers
for various differential equations. Equations within the realm of this package
include:

- Discrete equations (function maps, discrete stochastic (Gillespie/Markov)
  simulations)
- Ordinary differential equations (ODEs)
- Split and Partitioned ODEs (Symplectic integrators, IMEX Methods)
- Stochastic ordinary differential equations (SODEs or SDEs)
- Random differential equations (RODEs or RDEs)
- Differential algebraic equations (DAEs)
- Delay differential equations (DDEs)
- Mixed discrete and continuous equations (Hybrid Equations, Jump Diffusions)
- (Stochastic) partial differential equations ((S)PDEs) (with both finite
  difference and finite element methods)

The well-optimized DifferentialEquations solvers benchmark as the some of the fastest
implementations, using classic algorithms and ones from recent research which
routinely outperform the "standard" C/Fortran methods, and include algorithms
optimized for high-precision and HPC applications. At the same time, it wraps
the classic C/Fortran methods, making it easy to switch over to them whenever
necessary. Solving differential equations with different methods from
different languages and packages can be done by changing one line of code,
allowing for easy benchmarking to ensure you are using the fastest method possible.

DifferentialEquations.jl integrates with the Julia package sphere with:

- GPU accleration through CUDAnative.jl and CuArrays.jl
- Automated sparsity detection with [SparsityDetection.jl](@ref)
- Automatic Jacobian coloring with [SparseDiffTools.jl](@ref), allowing for fast solutions
  to problems with sparse or structured (Tridiagonal, Banded, BlockBanded, etc.) Jacobians
- Allowing the specification of linear solvers for maximal efficiency
- Progress meter integration with the Juno IDE for estimated time to solution
- Automatic plotting of time series and phase plots
- Built-in interpolations
- Wraps for common C/Fortran methods like Sundials and Hairer's radau
- Arbitrary precision with BigFloats and Arbfloats
- Arbitrary array types, allowing the definition of differential equations on
  matrices and distributed arrays
- Unit checked arithmetic with Unitful

Additionally, DifferentialEquations.jl comes with built-in analysis features, including:

- [Forward and adjoint local sensitivity analysis](@ref) for fast gradient computations
- [Optimization-based and Bayesian parameter estimation](@ref)
- Neural differential equations with [DiffEqFlux.jl](@ref)
  for efficient scientific machine learning (scientific ML) and scientific AI.
- [Automatic distributed, multithreaded, and GPU parallelism of ensemble trajectories](@ref)
- [Global sensitivity analysis](@ref)
- [Uncertainty quantification](@ref)

If you have any questions, or just want to chat about solvers/using the package,
please feel free to use the [Gitter channel](@ref).
For bug reports, feature requests, etc., please submit an issue. If you're
interested in contributing, please see the
[Developer Documentation](@ref).

## Supporting and Citing

The software in this ecosystem was developed as part of academic research.
If you would like to help support it, please star the repository as such
metrics may help us secure funding in the future. If you use JuliaDiffEq
software as part of your research, teaching, or other activities, we would
be grateful if you could cite our work.
[Please see our citation page for guidelines](@ref).

## Getting Started: Installation And First Steps

### Installing from Julia

To install the package, use the following command inside the Julia REPL:
```julia
using Pkg
Pkg.add("DifferentialEquations")
```

To load the package, use the command:

```julia
using DifferentialEquations
```

This will add solvers and dependencies
for all kinds of Differential Equations (e.g. ODEs or SDEs etc., see the Supported
Equations section below). If you are interested in only one type of equation
solvers of `DifferentialEquations.jl` or simply want a more lightweight
version, see the
[Low Dependency Usage](@ref)
page.

To understand the package in more detail, check out the following tutorials in
this manual. **It is highly recommended that new users start with the
[ODE tutorial](@ref ode_example)**. Example IJulia notebooks
[can also be found in DiffEqTutorials.jl](@ref).
If you find any example where there seems to be an error, please open an issue.

For the most up to date information on using the package, please join [the Gitter channel](@ref).

Using the bleeding edge for the latest features and development is only recommended
for power users. Information on how to get to the bleeding edge is found in the
[developer documentation](@ref).

### Installing from Python

Use of DifferentialEquations.jl from the Python programming language is available through the
[diffeqpy](@ref) module. To install diffeqpy, use pip:

```
pip install diffeqpy
```

Using diffeqpy requires that Julia is installed and in the path, along
with DifferentialEquations.jl and PyCall.jl. To install Julia,
download a generic binary from
[the JuliaLang site](@ref) and add it to your path.
To install Julia packages required for diffeqpy, open up Python
interpreter then run:

```pycon
>>> import diffeqpy
>>> diffeqpy.install()
```

and you're good! In addition, to improve the performance of your code it is
recommended that you use Numba to JIT compile your derivative functions. To
install Numba, use:

```
pip install numba
```

diffeqpy supports the majority of DifferentialEquations.jl with very similar
syntax, see [the diffeqpy README for more details](@ref).
One important point to note is that Numba is generally an order of magnitude slower
than Julia in terms of  the generated differential equation solver code, and thus it is
recommended to use `julia.Main.eval` for Julia-side derivative function implementations
for maximal efficiency. See [this blog post](@ref)
for more information.

### Installing from R

Use of DifferentialEquations.jl from the R programming language is available through the
[diffeqr](@ref) module.
[diffeqr is registered into CRAN](@ref).
Thus to add the package, use:

```R
install.packages("diffeqr")
```

To install the master branch of the package (for developers), use:

```R
devtools::install_github('JuliaDiffEq/diffeqr', build_vignettes=T)
```

You will need a working installation of Julia in your path. To install Julia, download a generic binary
from [the JuliaLang site](@ref) and add it to your path. The download and
installation of DifferentialEquations.jl will happen on the first invocation of `diffeqr::diffeq_setup()`.

Currently, use from R supported a subset of DifferentialEquations.jl which is documented
[through CRAN](@ref)

### IJulia Notebook Tutorials

You can access extra tutorials supplied in the
[DiffEqTutorials.jl repository](@ref)
via the commands:

```julia
using Pkg
pkg"add https://github.com/JuliaDiffEq/DiffEqTutorials.jl"
using DiffEqTutorials
DiffEqTutorials.open_notebooks()
```

Or you can view the webpages for the rendered tutorials at the links found in the repository.

### Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](@ref)

### Tutorials

The following tutorials will introduce you to the functionality of
DifferentialEquations.jl. More examples can be found by
[checking out the IJulia notebooks in the examples folder](@ref).

```@contents
Pages = [
    "tutorials/ode_example.md",
    "tutorials/sde_example.md",
    "tutorials/dde_example.md",
    "tutorials/dae_example.md",
    "tutorials/discrete_stochastic_example.md",
    "tutorials/jump_diffusion.md",
    "tutorials/bvp_example.md",
    "tutorials/additional.md"
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
    "basics/integrator.md",
    "basics/problem.md",
    "basics/faq.md",
    "basics/compatibility_chart.md"
    ]
Depth = 2
```


### Problem Types

These pages describe building the problem types to define differential equations
for the solvers, and the special features of the different solution types.

```@contents
Pages = [
  "types/discrete_types.md",
  "types/ode_types.md",
  "types/dynamical_types.md",
  "types/split_ode_types.md",
  "types/steady_state_types.md",
  "types/bvp_types.md",
  "types/sde_types.md",
  "types/rode_types.md",
  "types/dde_types.md",
  "types/dae_types.md",
  "types/jump_types.md",
]
Depth = 2
```

### Solver Algorithms

These pages describe the solvers and available algorithms in detail.

```@contents
Pages = [
  "solvers/discrete_solve.md",
  "solvers/ode_solve.md",
  "solvers/dynamical_solve.md",
  "solvers/split_ode_solve.md",
  "solvers/steady_state_solve.md",
  "solvers/bvp_solve.md",
  "solvers/jump_solve.md",
  "solvers/sde_solve.md",
  "solvers/rode_solve.md",
  "solvers/dde_solve.md",
  "solvers/dae_solve.md",
  "solvers/benchmarks.md"
]
Depth = 2
```

### Additional Features

These sections discuss extra performance enhancements, event handling, and other
in-depth features.

```@contents
Pages = [
    "features/performance_overloads.md",
    "features/diffeq_arrays.md",
    "features/diffeq_operator.md",
    "features/noise_process.md",
    "features/linear_nonlinear.md",
    "features/callback_functions.md",
    "features/callback_library.md",
    "features/ensemble.md",
    "features/io.md",
    "features/low_dep.md",
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
    "analysis/bifurcation.md",
    "analysis/sensitivity.md",
    "analysis/global_sensitivity.md",
    "analysis/uncertainty_quantification.md",
    "analysis/neural_networks.md",
    "analysis/dev_and_test.md"
]
Depth = 2
```

### Modeling Tools

While DifferentialEquations.jl can be used to directly build any differential
or difference equation (/ discrete stochastic) model, in many cases it can be
helpful to have a tailored-built API for making certain types of common models
easier. This is provided by the modeling functionality.

```@contents
Pages = [
    "models/multiscale.md",
    "models/physical.md",
    "models/financial.md",
    "models/biological.md",
    "models/external_modeling.md"
]
Depth = 2
```

### APIs

Some DifferentialEquations.jl packages provide documented APIs, these include:
```@contents
Pages = [
    "apis/diffeqbio.md"
]
Depth = 2
```



### Extra Details

These are just assorted extra explanations for the curious.

```@contents
Pages = [
    "extras/timestepping.md"
]
Depth = 2
```

## Acknowledgements

#### Core Contributors

JuliaDiffEq and DifferentialEquations.jl has been a collaborative effort by many
individuals. Significant contributions have been made by the following individuals:

- Chris Rackauckas (@ChrisRackauckas) (lead developer)
- Yingbo Ma (@YingboMa)
- David Widmann (@devmotion)
- Hendrik Ranocha (@ranocha)
- Ethan Levien (@elevien)
- Tom Short (@tshort)
- @dextorious
- Samuel Isaacson (@isaacsas)

#### Google Summer of Code Alumni

- Yingbo Ma (@YingboMa)
- Shivin Srivastava (@shivin9)
- Ayush Pandey (@Ayush-iitkgp)
- Xingjian Guo (@MSeeker1340)
- Shubham Maddhashiya (@sipah00)
- Vaibhav Kumar Dixit (@Vaibhavdixit02)
