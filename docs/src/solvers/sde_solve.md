# SDE Solvers

## Recommended Methods

For most problems where a good amount of accuracy is required and stiffness may
be an issue, the `SRIW1Optimized` algorithm should do well. If the problem
has additive noise, then `SRA1Optimized` will be the optimal algorithm. If
you simply need to quickly compute a large ensamble and don't need accuracy
(and don't have stiffness problems), then `EM` can do well.

## Special Keyword Arguments

* `discard_length` - Size at which to discard future information in adaptive. Default is 1e-15.
* `tableau`: The tableau for an `:SRA` or `:SRI` algorithm. Defaults to SRIW1 or SRA1.
* `adaptivealg`: The adaptive timestepping algorithm. Default is `:RSwm3`.

## Implemented Solvers

### StochasticDiffEq.jl

- `EM`- The Euler-Maruyama method.
- `RKMil` - An explicit Runge-Kutta discretization of the strong Order 1.0 Milstein method.
- `SRA` - The strong Order 2.0 methods for additive SDEs due to Rossler. Not yet implemented.
  Default tableau is for SRA1.
- `SRI` - The strong Order 1.5 methods for diagonal/scalar SDEs due to Rossler.
  Default tableau is for SRIW1.
- `SRIW1` - An optimized version of SRIW1. Strong Order 1.5.
- `SRA1` - An optimized version of SRIA1. Strong Order 2.0.
