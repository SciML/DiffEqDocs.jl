# FEM Poisson Solvers

## Recommended Methods

The only available method is `FEMDiffEqPoisson`. This method uses a chosen linear
solver from IterativeSolvers.jl for a linear problem or a nonlinear solver
from NLsolve.jl for a nonlinear problem.

# Full List of Methods

## FiniteElementDiffEq.jl

* Factorizations (`:LU`, `:Cholesky`, `:QR`, `:SVD`)
* Conjugate-Gradient (`:CG`)
* `:GMRES`

Example:

```julia
sol = solve(prob,FEMDiffEqPoisson(),solver=:CG)
```
