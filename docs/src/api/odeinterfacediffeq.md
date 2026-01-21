# [ODEInterfaceDiffEq.jl](@id odeinterfacediffeq_api)

ODEInterfaceDiffEq.jl provides a wrapper for classic Fortran ODE solvers from
the ODEInterface.jl package into the SciML interface. These include the Hairer
and Wanner solvers.

Note that this package is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use ODEInterfaceDiffEq.jl:

```julia
using Pkg
Pkg.add("ODEInterfaceDiffEq")
import ODEInterfaceDiffEq
```

These methods can be used independently of the rest of DifferentialEquations.jl.

## ODE Solver APIs

```@docs
ODEInterfaceDiffEq.dopri5
ODEInterfaceDiffEq.dop853
ODEInterfaceDiffEq.odex
ODEInterfaceDiffEq.seulex
ODEInterfaceDiffEq.radau
ODEInterfaceDiffEq.radau5
ODEInterfaceDiffEq.rodas
ODEInterfaceDiffEq.ddeabm
ODEInterfaceDiffEq.ddebdf
```

## BVP Solver APIs

```@docs
ODEInterfaceDiffEq.BVPM2
ODEInterfaceDiffEq.BVPSOL
ODEInterfaceDiffEq.COLNEW
```
