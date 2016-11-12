# Stokes Types

## Problem Type

`StokesProblem`
Defines the solution to a stationary Stokes problem:
```math
```
### Constructors
`StokesProblem(f₁,f₂,g,uanalytic,vanalytic,panalytic)`
`StokesProblem(f₁,f₂,g,ugD,vgD)`
### Fields
* `f₁::Function`
* `f₂::Function`
* `g::Function`
* `ugD::Function`
* `vgD::Function`
* `uanalytic::Function`
* `vanalytic::Function`
* `panalytic::Function`
* `trueknown::Bool`

## Example Problems

Examples problems can be found in [src/premades/premade_problems.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl/blob/master/src/premades/premade_problems.jl)
