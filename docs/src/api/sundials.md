# [Sundials.jl](@id sundials)

This is a wrapper package for importing solvers from Sundials into the SciML interface.
Note that these solvers do not come by default, and thus one needs to install
the package before using these solvers:

```julia
]add Sundials
using Sundials
```

These methods can be used independently of the rest of DifferentialEquations.jl.

## ODE Solver APIs

```@docs
CVODE_Adams
CVODE_BDF
ARKODE
```

## DAE Solver APIs

```@docs
IDA
```
