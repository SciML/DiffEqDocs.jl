# [BVP Solvers](@id bvp_solve)

```julia
solve(prob::BVProblem, alg, dt; kwargs)
solve(prob::TwoPointBVProblem, alg, dt; kwargs)
solve(prob::SecondOrderBVProblem, alg, dt; kwargs)
```

Solves the BVP defined by `prob` using the algorithm `alg`. All algorithms except `Shooting` and `MultipleShooting` methods should specify a `dt` which is the step size for the discretized mesh.

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package | Methods | Good for |
|---|---|---|
| `BoundaryValueDiffEqShooting` | `Shooting`, `MultipleShooting` | Fastest when not too stiff. |
| `BoundaryValueDiffEqMIRK` | MIRK2/3/4/5/6 (sparse Jacobians) | Default BVP workhorse; robust on standard two-point BVPs. |
| `BoundaryValueDiffEqFIRK` | RadauIIa1/2/3/5/7, LobattoIIIa/b/c | Stiff or high-precision BVPs (fully-implicit RK). |
| `BoundaryValueDiffEqAscher` | Ascher1/2/3/4/5/6/7 | Index-1 DAEs and mixed-order / stiff BVPs. |
| `BoundaryValueDiffEqMIRKN` | MIRK-N methods | Second-order BVPs (e.g. structural / mechanics). |
| `SimpleBoundaryValueDiffEq` | SimpleMIRK*, SimpleShooting | Lightweight BVP shooting; minimal dependencies. |
| `ODEInterface` | `BVPM2`, `BVPSOL` | Fortran BVP wrappers from ODEInterface.jl. |


## Recommended Methods

The `MIRK` methods are recommended in most scenarios given their improved stability properties over the other methods. They have adaptivty and sparsity handling which allows for them to handle large-scale and difficult problems. However, they are not compatible with callbacks / event handling (i.e. discontinuities), and in such cases [Shooting methods](https://en.wikipedia.org/wiki/Shooting_method) are required. There are single shooting methods and multiple shooting methods available in BoundaryValueDiffEq.jl. Shooting methods should be used with an appropriate ODE solver such as `Shooting(Tsit5())` or `MultipleShooting(5, FBDF())`. Additionally,
in many cases, single shooting method `Shooting` may be faster than collocation methods if it converges, though it is a lot less numerically robust. Multiple shooting method `MultipleShooting` is more stable and robust than single shooting method `Shooting`.

## Full List of Methods

### BoundaryValueDiffEq.jl

!!! note "v8: BoundaryValueDiffEq must be loaded explicitly"

    Under DifferentialEquations.jl v8 the `using DifferentialEquations`
    umbrella only re-exports `OrdinaryDiffEq`. The BVP solvers below come from
    `BoundaryValueDiffEq.jl` (and its sublibraries). Either load the umbrella

    ```julia
    using BoundaryValueDiffEq         # MIRK*, Shooting, MultipleShooting,
                                      # RadauIIa*, LobattoIIIa/b/c*, Ascher*,
                                      # MIRKN*, BVPJacobianAlgorithm
    ```

    or, for tighter compile times, just the relevant sublibrary:

    | Method family                                              | Sublibrary                          |
    |------------------------------------------------------------|-------------------------------------|
    | `Shooting`, `MultipleShooting`                             | `BoundaryValueDiffEqShooting`       |
    | `MIRK2`–`MIRK6`                                            | `BoundaryValueDiffEqMIRK`           |
    | `MIRKN4`, `MIRKN6`                                         | `BoundaryValueDiffEqMIRKN`          |
    | `LobattoIIIa*`, `LobattoIIIb*`, `LobattoIIIc*`, `RadauIIa*` | `BoundaryValueDiffEqFIRK`           |
    | `Ascher1`–`Ascher7` (BVDAE)                                | `BoundaryValueDiffEqAscher`         |
    | `BVPM2`, `BVPSOL`, `COLNEW` (Fortran wrappers)             | `BoundaryValueDiffEq` + `ODEInterface` |

    Shooting / `MultipleShooting` also need an OrdinaryDiffEq inner solver, e.g.
    `Shooting(Tsit5())` requires `using OrdinaryDiffEqTsit5: Tsit5`.

#### Shooting Methods

  - `BoundaryValueDiffEqShooting.Shooting(odealg())` - A wrapper over initial value problem solvers, it reduces BVP to an initial value problem and solves the IVP.
  - `BoundaryValueDiffEqShooting.MultipleShooting(N, odealg())` - A wrapper over initial value problem solvers, it reduces BVP to `N` initial value problems and solves these IVPs. Multiple Shooting usually maintains more numerical stability than Single Shooting.

#### MIRK(Monotonic Implicit Runge-Kutta) Methods

All `MIRK` methods have defect control adaptivity by default which adapts the mesh (`dt`) automatically. This can be turned
off via the keyword argument `adaptive = false`.

  - `BoundaryValueDiffEqMIRK.MIRK2` - A 2nd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `BoundaryValueDiffEqMIRK.MIRK3` - A 3rd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `BoundaryValueDiffEqMIRK.MIRK4` - A 4th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `BoundaryValueDiffEqMIRK.MIRK5` - A 5th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `BoundaryValueDiffEqMIRK.MIRK6` - A 6th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.

#### FIRK(Fully Implicit Runge-Kutta) methods

Similar to `MIRK` methods, fully implicit Runge-Kutta methods construct nonlinear problems from the collocation equations of a BVP and solve such nonlinear systems to obtain numerical solutions of BVP. When solving large boundary value problems, choose a nested NonlinearSolve.jl solver by setting `nested_nlsolve=true` in FIRK solvers can achieve better performance.

  - `BoundaryValueDiffEqFIRK.LobattoIIIa2` - A 2nd stage LobattoIIIa collocation method.

  - `BoundaryValueDiffEqFIRK.LobattoIIIa3` - A 3rd stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIa4` - A 4th stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIa5` - A 5th stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIb2` - A 2nd stage LobattoIIIa collocation method, doesn't support defect control adaptivity.
  - `BoundaryValueDiffEqFIRK.LobattoIIIb3` - A 3rd stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIb4` - A 4th stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIb5` - A 5th stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIc2` - A 2nd stage LobattoIIIa collocation method, doesn't support defect control adaptivity.
  - `BoundaryValueDiffEqFIRK.LobattoIIIc3` - A 3rd stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIc4` - A 4th stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.LobattoIIIc5` - A 5th stage LobattoIIIa collocation method.
  - `BoundaryValueDiffEqFIRK.RadauIIa1` - A 1st stage Radau collocation method, doesn't support defect control adaptivity.
  - `BoundaryValueDiffEqFIRK.RadauIIa2` - A 2nd stage Radau collocation method.
  - `BoundaryValueDiffEqFIRK.RadauIIa3` - A 3rd stage Radau collocation method.
  - `BoundaryValueDiffEqFIRK.RadauIIa5` - A 5th stage Radau collocation method.
  - `BoundaryValueDiffEqFIRK.RadauIIa7` - A 7th stage Radau collocation method.

#### Gauss Legendre collocation methods

The `Ascher` collocation methods are similar with `MIRK` and `FIRK` methods but have extension for BVDAE prblem solving, the error control is based on instead of defect control adaptivity.

  - `BoundaryValueDiffEqAscher.Ascher1` - A 1st stage Gauss Legendre collocation method with Ascher's error control adaptivity.
  - `BoundaryValueDiffEqAscher.Ascher2` - A 2nd stage Gauss Legendre collocation method with Ascher's error control adaptivity.
  - `BoundaryValueDiffEqAscher.Ascher3` - A 3rd stage Gauss Legendre collocation method with Ascher's error control adaptivity.
  - `BoundaryValueDiffEqAscher.Ascher4` - A 4th stage Gauss Legendre collocation method with Ascher's error control adaptivity.
  - `BoundaryValueDiffEqAscher.Ascher5` - A 5th stage Gauss Legendre collocation method with Ascher's error control adaptivity.
  - `BoundaryValueDiffEqAscher.Ascher6` - A 6th stage Gauss Legendre collocation method with Ascher's error control adaptivity.
  - `BoundaryValueDiffEqAscher.Ascher7` - A 7th stage Gauss Legendre collocation method with Ascher's error control adaptivity.

#### MIRKN(Monotonic Implicit Runge-Kutta-Nystöm) methods

  - `BoundaryValueDiffEqMIRKN.MIRKN4` - A 4th order collocation method using an implicit Runge-Kutta-Nyström tableau without defect control adaptivity.
  - `BoundaryValueDiffEqMIRKN.MIRKN6` - A 6th order collocation method using an implicit Runge-Kutta-Nyström tableau without defect control adaptivity.

### SimpleBoundaryValueDiffEq.jl

  - `SimpleBoundaryValueDiffEq.SimpleMIRK4` - A simplified 4th order collocation method using an implicit Runge-Kutta tableau.
  - `SimpleBoundaryValueDiffEq.SimpleMIRK5` - A simplified 5th order collocation method using an implicit Runge-Kutta tableau.
  - `SimpleBoundaryValueDiffEq.SimpleMIRK6` - A simplified 6th order collocation method using an implicit Runge-Kutta tableau.
  - `SimpleBoundaryValueDiffEq.SimpleShooting` - A simplified single Shooting method.

### ODEInterface.jl

ODEInterface.jl can be used seamlessly with BoundaryValueDiffEq.jl, after we define our model using `BVProblem` or `TwoPointBVProblem`, we can directly call the solvers from ODEInterface.jl.

  - `ODEInterface.BVPM2` - FORTRAN code for solving two-point boundary value problems. `BVPM2` is only compatible with `TwoPointBVProblem`.
  - `ODEInterface.BVPSOL` - FORTRAN77 code which solves highly nonlinear two point boundary value problems using a local linear solver (condensing algorithm) or a global sparse linear solver for the solution of the arising linear subproblems, by Peter Deuflhard, Georg Bader, Lutz Weimann. `BVPSOL` should be used with `TwoPointBVProblem` and initial guess.
  - `ODEInterface.COLNEW` - A Fortran77 code solves a multi-points boundary value problems for a mixed order system of ODEs by Uri Ascher and Georg Bader. It incorporates a new basis representation replacing b-splines, and improvements for the linear and nonlinear algebraic equation solvers. `COLNEW` support `TwoPointBVProblem` by default. To solve multi-points BVP using `COLNEW`, special form of multi-points boundary conditions should be provided by `COLNEW(bc_func, dbc_func, zeta)` where `bc_func(i, z, res)` is the multi-points boundary conditions, `dbc_func(i, z, dbc)` is the i-th row of jacobian of boundary conditions.
