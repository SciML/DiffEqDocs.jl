# [SteadyStateDiffEq.jl](@id steadystatediffeq_api)

SteadyStateDiffEq.jl is the native Julia package for solving steady state problems
within the SciML ecosystem. It provides methods for finding equilibrium solutions
of differential equations.

## Installation

SteadyStateDiffEq.jl is included with DifferentialEquations.jl. To use it standalone:

```julia
using Pkg
Pkg.add("SteadyStateDiffEq")
import SteadyStateDiffEq
```

## Steady State Solver APIs

```@docs
SteadyStateDiffEq.DynamicSS
SteadyStateDiffEq.SSRootfind
```
