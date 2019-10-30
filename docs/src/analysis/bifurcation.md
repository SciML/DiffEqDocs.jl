# Bifurcation Analysis

Bifurcation analysis on DifferentialEquations.jl types can be performed by:

- Bifurcations.jl
- PseudoArcLengthContinuations.jl
- PyDSTool.jl
- DiffEqBiological.jl

If your system is a chemical reaction system, see the documentation at 
[DiffEqBiological.jl](https://github.com/JuliaDiffEq/DiffEqBiological.jl#making-bifurcation-diagram)
for quickly generating bifurcation plots. Bifurcations.jl can directly generate a `BifurcationProblem`
from an `ODEProblem`. PyDSTool.jl is no longer recommended.
