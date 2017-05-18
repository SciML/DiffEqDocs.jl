# FEM Heat Solvers

## Recommended Methods

For nonstiff problems it's recommended you use `FEMDiffEqHeatEuler`, while for
stiff problems it's recommended that you use `FEMDiffEqHeatSemiImplicitCrankNicholson`.

## Avaliable Methods

* `FEMDiffEqHeatEuler`
* `FEMDiffEqHeatImplicitEuler`
* `FEMDiffEqHeatCrankNicholson`
* `FEMDiffEqHeatSemiImplicitEuler`
* `FEMDiffEqHeatSemiImplicitCrankNicholson`

Additionally, for linear solves, one can choose the method by which the linear solve
takes place via the `method` keyword argument.

* Factorizations (`:LU`, `:Cholesky`, `:QR`, `:SVD`)
* Conjugate-Gradient (`:CG`)
* `:GMRES`

Example:

```julia
sol = solve(prob,FEMDiffEqHeatCrankNicholson(),solver=:CG)
```
