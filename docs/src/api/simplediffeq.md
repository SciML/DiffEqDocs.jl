# [SimpleDiffEq.jl](@id simplediffeq_api)

SimpleDiffEq.jl provides simplified implementations of a few ODE and SDE solvers.
They are primarily designed for experimentation and offer shorter compile times.
They have limitations compared to OrdinaryDiffEq.jl and StochasticDiffEq.jl and
are generally not faster.

Note that this package is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use SimpleDiffEq.jl:

```julia
using Pkg
Pkg.add("SimpleDiffEq")
import SimpleDiffEq
```

## ODE Solver APIs

```@docs
SimpleDiffEq.SimpleTsit5
SimpleDiffEq.SimpleATsit5
SimpleDiffEq.GPUSimpleATsit5
SimpleDiffEq.SimpleEuler
SimpleDiffEq.SimpleRK4
SimpleDiffEq.GPUVern7
SimpleDiffEq.GPUVern9
```

## SDE Solver APIs

```@docs
SimpleDiffEq.SimpleEM
```
