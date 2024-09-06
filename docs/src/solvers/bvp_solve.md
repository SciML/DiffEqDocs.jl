# BVP Solvers

`solve(prob::BVProblem,alg,dt=0.0;kwargs)`

Solves the BVP defined by `prob` using the algorithm `alg`. All algorithms except `Shooting` methods should specify a `dt` which is the step size for the discretized mesh.

## Recommended Methods

The `MIRK` methods are recommended in most scenarios given their improved stability properties over the other methods. They have adaptivty and sparsity handling which allows for them to handle large-scale and difficult problems. However, they are not compatible with callbacks / event handling (i.e. discontinuities), and in such cases [Shooting methods](https://en.wikipedia.org/wiki/Shooting_method) are required. There are single shooting methods and multiple shooting methods available in BoundaryValueDiffEq.jl. Shooting methods should be used with an appropriate ODE solver such as `Shooting(Tsit5())` or `MultipleShooting(5, FBDF())`. Additionally,
in many cases, single shooting method `Shooting` may be faster than collocation methods if it converges, though it is a lot less numerically robust. Multiple shooting method `MultipleShooting` is more stable and robust than single shooting method `Shooting`.

## Full List of Methods

### BoundaryValueDiffEq.jl

#### Shooting Methods

  - `Shooting(odealg())` - A wrapper over initial value problem solvers, it reduces BVP to an initial value problem and solves the IVP.
  - `MultipleShooting(N, odealg())` - A wrapper over initial value problem solvers, it reduces BVP to `N` initial value problems and solves these IVPs. Multiple Shooting usually maintains more numerical stability than Single Shooting.

#### MIRK(Monotonic Implicit Runge-Kutta) Methods

All `MIRK` methods have defect control adaptivity by default which adapts the mesh (`dt`) automatically. This can be turned
off via the keyword argument `adaptive = false`.

  - `MIRK2` - A 2nd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK3` - A 3rd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK4` - A 4th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK5` - A 5th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK6` - A 6th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.

#### FIRK(Fully Implicit Runge-Kutta) methods

Similar to `MIRK` methods, fully implicit Runge-Kutta methods construct nonlinear problems from the collocation equations of a BVP and solve such nonlinear systems to obtain numerical solutions of BVP. When solving large boundary value problems, choose a nested NonlinearSolve.jl solver by setting `nested_nlsolve=true` in FIRK solvers can achieve better performance.

  - `LobattoIIIa2` - A 2nd stage LobattoIIIa collocation method.

  - `LobattoIIIa3` - A 3rd stage LobattoIIIa collocation method.
  - `LobattoIIIa4` - A 4th stage LobattoIIIa collocation method.
  - `LobattoIIIa5` - A 5th stage LobattoIIIa collocation method.
  - `LobattoIIIb2` - A 2nd stage LobattoIIIa collocation method, doesn't support defect control adaptivity.
  - `LobattoIIIb3` - A 3rd stage LobattoIIIa collocation method.
  - `LobattoIIIb4` - A 4th stage LobattoIIIa collocation method.
  - `LobattoIIIb5` - A 5th stage LobattoIIIa collocation method.
  - `LobattoIIIc2` - A 2nd stage LobattoIIIa collocation method, doesn't support defect control adaptivity.
  - `LobattoIIIc3` - A 3rd stage LobattoIIIa collocation method.
  - `LobattoIIIc4` - A 4th stage LobattoIIIa collocation method.
  - `LobattoIIIc5` - A 5th stage LobattoIIIa collocation method.
  - `RadauIIa1` - A 1st stage Radau collocation method, doesn't support defect control adaptivity.
  - `RadauIIa2` - A 2nd stage Radau collocation method.
  - `RadauIIa3` - A 3rd stage Radau collocation method.
  - `RadauIIa5` - A 5th stage Radau collocation method.
  - `RadauIIa7` - A 7th stage Radau collocation method.

### ODEInterface.jl

ODEInterface.jl can be used seamlessly with BoundaryValueDiffEq.jl, after we define our model using `BVProblem` or `TwoPointBVProblem`, we can directly call the solvers from ODEInterface.jl.

  - `BVPM2` - FORTRAN code for solving two-point boundary value problems. `BVPM2` is only compatible with `TwoPointBVProblem`.
  - `BVPSOL` - FORTRAN77 code which solves highly nonlinear two point boundary value problems using a local linear solver (condensing algorithm) or a global sparse linear solver for the solution of the arising linear subproblems, by Peter Deuflhard, Georg Bader, Lutz Weimann. `BVPSOL` should be used with `TwoPointBVProblem` and initial guess.
  - `COLNEW` - A Fortran77 code solves a multi-points boundary value problems for a mixed order system of ODEs by Uri Ascher and Georg Bader. It incorporates a new basis representation replacing b-splines, and improvements for the linear and nonlinear algebraic equation solvers. `COLNEW` support `TwoPointBVProblem` by default. To solve multi-points BVP using `COLNEW`, special form of multi-points boundary conditions should be provided by `COLNEW(bc_func, dbc_func, zeta)` where `bc_func(i, z, res)` is the multi-points boundary conditions, `dbc_func(i, z, dbc)` is the i-th row of jacobian of boundary conditions.
