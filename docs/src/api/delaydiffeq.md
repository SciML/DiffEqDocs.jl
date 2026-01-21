# [DelayDiffEq.jl](@id delaydiffeq_api)

DelayDiffEq.jl is the native Julia package for solving delay differential equations
(DDEs) within the SciML ecosystem. It provides methods for constant and state-dependent
delays using the method of steps.

## Installation

DelayDiffEq.jl is included with DifferentialEquations.jl. To use it standalone:

```julia
using Pkg
Pkg.add("DelayDiffEq")
import DelayDiffEq
```

## DDE Solver API

```@docs
DelayDiffEq.MethodOfSteps
```
