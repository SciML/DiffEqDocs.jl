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
BoundaryValueDiffEq.Shooting
BoundaryValueDiffEq.MultipleShooting
```

## MIRK Method APIs

```@docs
BoundaryValueDiffEq.MIRK2
BoundaryValueDiffEq.MIRK3
BoundaryValueDiffEq.MIRK4
BoundaryValueDiffEq.MIRK5
BoundaryValueDiffEq.MIRK6
```

## FIRK Method APIs (Lobatto)

```@docs
BoundaryValueDiffEq.LobattoIIIa2
BoundaryValueDiffEq.LobattoIIIa3
BoundaryValueDiffEq.LobattoIIIa4
BoundaryValueDiffEq.LobattoIIIa5
BoundaryValueDiffEq.LobattoIIIb2
BoundaryValueDiffEq.LobattoIIIb3
BoundaryValueDiffEq.LobattoIIIb4
BoundaryValueDiffEq.LobattoIIIb5
BoundaryValueDiffEq.LobattoIIIc2
BoundaryValueDiffEq.LobattoIIIc3
BoundaryValueDiffEq.LobattoIIIc4
BoundaryValueDiffEq.LobattoIIIc5
```

## FIRK Method APIs (Radau)

```@docs
BoundaryValueDiffEq.RadauIIa1
BoundaryValueDiffEq.RadauIIa2
BoundaryValueDiffEq.RadauIIa3
BoundaryValueDiffEq.RadauIIa5
BoundaryValueDiffEq.RadauIIa7
```

## Ascher Collocation Method APIs

```@docs
BoundaryValueDiffEq.Ascher1
BoundaryValueDiffEq.Ascher2
BoundaryValueDiffEq.Ascher3
BoundaryValueDiffEq.Ascher4
BoundaryValueDiffEq.Ascher5
BoundaryValueDiffEq.Ascher6
BoundaryValueDiffEq.Ascher7
```

## MIRKN Method APIs

```@docs
BoundaryValueDiffEq.MIRKN4
BoundaryValueDiffEq.MIRKN6
```
