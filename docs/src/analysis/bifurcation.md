# Bifurcation Analysis

Bifurcation analysis on DifferentialEquations.jl types can be performed by:

- [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl) (currently the most comprehensive and activately maintained package)
- [Bifurcations.jl](https://github.com/tkf/Bifurcations.jl)
- PyDSTool.jl (no longer recommended)

BifurcationKit has integration with `ODEProblem` for some functionality (like computing periodic orbits via shooting).
If `oprob` is an `ODEProblem`, one can also just pass `oprob.f.f` and `oprob.f.jac` to BifurcationKit methods as needed.
Bifurcations.jl can directly generate a `BifurcationProblem` from an
`ODEProblem`. 
