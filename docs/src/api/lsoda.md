# [LSODA.jl](@id lsoda_api)

LSODA.jl is a wrapper package that provides access to the LSODA algorithm within
the SciML interface. LSODA is a well-known method which uses automatic switching
to solve both stiff and non-stiff equations.

Note that this package is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use LSODA.jl:

```julia
using Pkg
Pkg.add("LSODA")
import LSODA
```

These methods can be used independently of the rest of DifferentialEquations.jl.

## ODE Solver APIs

```@docs
LSODA.lsoda
```
