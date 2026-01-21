# [DASSL.jl](@id dassl_api)

DASSL.jl is a native Julia implementation of the DASSL algorithm for solving
differential-algebraic equations (DAEs) within the SciML interface.

Note that this package is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use DASSL.jl:

```julia
using Pkg
Pkg.add("DASSL")
import DASSL
```

These methods can be used independently of the rest of DifferentialEquations.jl.

## DAE Solver APIs

```@docs
DASSL.dassl
```
