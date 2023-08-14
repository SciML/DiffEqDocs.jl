# BVP Solvers

## Recommended Methods

### Shooting method

`Shooting` method paired with an appropriate integrator for the IVP, such as
`Shooting(Tsit5())`, is a flexible and efficient option to solve boundary value problems.
This allows one to combine callbacks/event handling with the BVP solver, and the high-order
interpolations can be used to define complex boundary conditions. However,
`Shooting` methods can in some cases be prone to sensitivity of the boundary
condition and its numerical results cannot always be trusted.

### MIRK method

`MIRK` methods utilize a sparse Jacobian to greatly improve efficiency. Besides, with defect control adaptivity turned on by default(to turn adaptivity, we can simply `solve(prob, MIRK4(), adaptivity=false)`), `MIRK` methods can automatically optimize the mesh, reduce computing cost and improve solution accuracy.

## Full List of Methods

### BoundaryValueDiffEq.jl

  - `Shooting` - A wrapper over initial value problem solvers, it reduces bvp to an initial value problem and solves the ivp.
  - `MIRK3` - A 3rd order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian
  - `MIRK4` - A 4th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK5` - A 5th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
  - `MIRK6` - A 6th order collocation method using an implicit Runge-Kutta tableau with a sparse Jacobian.
