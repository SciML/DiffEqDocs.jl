# FEM Poisson Solvers

## Recommended Methods

The only available method is `FEMDiffEqPoisson`. This method uses a chosen linear
solver from IterativeSolvers.jl for a linear problem or a nonlinear solver
from NLsolve.jl for a nonlinear problem. Additionally, the keyword `method` can
be used to specify the

### List of Methods

* Factorizations (`:LU`, `:Cholesky`, `:QR`, `:SVD`)
* Conjugate-Gradient (`:CG`)
* `:GMRES`

Example:

```julia
sol = solve(prob,FEMDiffEqPoisson(),solver=:CG)
```
