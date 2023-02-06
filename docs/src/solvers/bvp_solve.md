# BVP Solvers

## Recommended Methods

`GeneralMIRK4` is a well-rounded method when the problem is not too large.
It uses highly stable trust region methods to solve a 4th order fully
implicit Runge-Kutta scheme. As an alternative on general `BVProblem`s, the
`Shooting` method paired with an appropriate integrator for the IVP, such as
`Shooting(Tsit5())`, is a flexible and efficient option. This allows one
to combine callbacks/event handling with the BVP solver, and the high order
interpolations can be used to define complex boundary conditions. However,
`Shooting` methods can in some cases be prone to sensitivity of the boundary
condition.

When the problem is a large two-point boundary value problem that is sensitive
to the boundary conditions, `MIRK4` utilizes a sparse Jacobian to greatly
improve the efficiency.

## Full List of Methods

### BoundaryValueDiffEq.jl

  - `Shooting` - A wrapper over initial value problem solvers.
  - `GeneralMIRK4` - A 4th order collocation method using an implicit Runge-Kutta
    tableau, solved using a trust region dogleg method from NLsolve.jl.
  - `MIRK4` - A 4th order collocation method using an implicit Runge-Kutta tableau
    with a sparse Jacobian. Compatible only with two-point boundary value problems.
