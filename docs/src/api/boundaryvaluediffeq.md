# [BoundaryValueDiffEq.jl](@id boundaryvaluediffeq_api)

BoundaryValueDiffEq.jl is the native Julia package for solving boundary value problems
(BVPs) within the SciML ecosystem. It provides shooting methods, MIRK collocation methods,
FIRK methods, and more.

## Installation

BoundaryValueDiffEq.jl is included with DifferentialEquations.jl. To use it standalone:

```julia
using Pkg
Pkg.add("BoundaryValueDiffEq")
import BoundaryValueDiffEq
```

## Shooting Method APIs

```@docs
BoundaryValueDiffEqShooting.Shooting
BoundaryValueDiffEqShooting.MultipleShooting
```

## MIRK Method APIs

```@docs
BoundaryValueDiffEqMIRK.MIRK2
BoundaryValueDiffEqMIRK.MIRK3
BoundaryValueDiffEqMIRK.MIRK4
BoundaryValueDiffEqMIRK.MIRK5
BoundaryValueDiffEqMIRK.MIRK6
```

## FIRK Method APIs (Lobatto)

```@docs
BoundaryValueDiffEqFIRK.LobattoIIIa2
BoundaryValueDiffEqFIRK.LobattoIIIa3
BoundaryValueDiffEqFIRK.LobattoIIIa4
BoundaryValueDiffEqFIRK.LobattoIIIa5
BoundaryValueDiffEqFIRK.LobattoIIIb2
BoundaryValueDiffEqFIRK.LobattoIIIb3
BoundaryValueDiffEqFIRK.LobattoIIIb4
BoundaryValueDiffEqFIRK.LobattoIIIb5
BoundaryValueDiffEqFIRK.LobattoIIIc2
BoundaryValueDiffEqFIRK.LobattoIIIc3
BoundaryValueDiffEqFIRK.LobattoIIIc4
BoundaryValueDiffEqFIRK.LobattoIIIc5
```

## FIRK Method APIs (Radau)

```@docs
BoundaryValueDiffEqFIRK.RadauIIa1
BoundaryValueDiffEqFIRK.RadauIIa2
BoundaryValueDiffEqFIRK.RadauIIa3
BoundaryValueDiffEqFIRK.RadauIIa5
BoundaryValueDiffEqFIRK.RadauIIa7
```

## Ascher Collocation Method APIs

```@docs
BoundaryValueDiffEqAscher.Ascher1
BoundaryValueDiffEqAscher.Ascher2
BoundaryValueDiffEqAscher.Ascher3
BoundaryValueDiffEqAscher.Ascher4
BoundaryValueDiffEqAscher.Ascher5
BoundaryValueDiffEqAscher.Ascher6
BoundaryValueDiffEqAscher.Ascher7
```

## MIRKN Method APIs

```@docs
BoundaryValueDiffEqMIRKN.MIRKN4
BoundaryValueDiffEqMIRKN.MIRKN6
```
