# BVP Solvers

## Recommended Methods

The `MIRK` methods are recommended in most scenarios given their improved stability properties
over the other methods. They have adaptivty and sparsity handling which allows for them to
handle large-scale and difficult problems. However, they are not compatible with callbacks / event handling
(i.e. discontinuities), and in such cases `Shooting` is required, where `Shooting` should be
used with an appropriate ODE solver such as `Shooting(Tsit5())` or `Shooting(FBDF())`. Additionally,
in many cases `Shooting` may be faster if it converges, though it is a lot less numerically robust.

## Full List of Methods

### BoundaryValueDiffEq.jl

#### Shooting Methods

  - `Shooting` - A wrapper over initial value problem solvers, it reduces bvp to an initial value problem and solves the ivp.

#### MIRK Collocation Methods

All `MIRK` methods have defect control adaptivity by default which adapts the mesh (`dt`) automatically. This can be turned
off via the keyword argument `adaptive = false`.

  - `MIRK2` - A 2nd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK3` - A 3rd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK4` - A 4th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK5` - A 5th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK6` - A 6th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
