# DASKR.jl

This is a wrapper package for importing solvers from DASKR into the SciML interface.
DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
]add DASKR
using DASKR
```

These methods can be used independently of the rest of DifferentialEquations.jl.

## DAE Solver APIs

```@docs
daskr
```
