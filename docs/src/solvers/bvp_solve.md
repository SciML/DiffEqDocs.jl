# BVP Solvers

## Recommended Methods

The `MIRK` methods are recommended in most scenarios given their improved stability properties over the other methods. They have adaptivty and sparsity handling which allows for them to handle large-scale and difficult problems. However, they are not compatible with callbacks / event handling (i.e. discontinuities), and in such cases [Shooting methods](https://en.wikipedia.org/wiki/Shooting_method) is required. There are single shooting methods and multiple shooting methods available in BoundaryValueDiffEq.jl. Shooting methods should be used with an appropriate ODE solver such as `Shooting(Tsit5())` or `MultipleShooting(5, FBDF())`. Additionally,
in many cases, single shooting method `Shooting` may be faster if it converges, though it is a lot less numerically robust. Multiple shooting method `MultipleShooting` is more stable and robust than single shooting method `Shooting`.

## Full List of Methods

### BoundaryValueDiffEq.jl

#### Shooting Methods

  - `Shooting(odealg())` - A wrapper over initial value problem solvers, it reduces BVP to an initial value problem and solves the IVP.
  - `MultipleShooting(N, odealg())` - A wrapper over initial value problem solvers, it reduces BVP to `N` initial value problems and solves these IVPs. Multiple Shooting usually maintains more numerical stability than Single Shooting.

#### MIRK Collocation Methods

All `MIRK` methods have defect control adaptivity by default which adapts the mesh (`dt`) automatically. This can be turned
off via the keyword argument `adaptive = false`.

  - `MIRK2` - A 2nd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK3` - A 3rd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK4` - A 4th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK5` - A 5th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK6` - A 6th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.

### ODEInterface.jl

ODEInterface.jl can be used seamlessly with BoundaryValueDiffEq.jl, after we define our model using `BVProblem` or `TwoPointBVProblem`, we can directly call the solvers from ODEInterface.jl.

  - `BVPM2` - FORTRAN code for solving two-point boundary value problems. `BVPM2` is only compatiable with `TwoPointBVProblem`.
  - `BVPSOL` - FORTRAN77 code which solves highly nonlinear two point boundary value problems using a local linear solver (condensing algorithm) or a global sparse linear solver for the solution of the arising linear subproblems, by Peter Deuflhard, Georg Bader, Lutz Weimann. `BVPSOL` should be used with `TwoPointBVProblem` and initial guess.