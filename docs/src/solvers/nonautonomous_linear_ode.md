# Non-autonomous Linear ODE / Lie Group ODE Solvers

Non-autonomous linear ODE solvers focus on equations in the general form of

```math
u^\prime = A(u,p,t)u
```

and utilize the Lie group structure of the solution to accelerate the numerical
methods and capture certain properties in the solution process. One common simplification
is for solvers to require state-independent operators, which implies the form:

```math
u^\prime = A(t)u
```

Another type of solvers are needed when the operators are state-dependent, i.e.

```math
u^\prime = A(u)u
```

Others specifically require linearity, i.e.

```math
u^\prime = Au
```

where ``A`` is a constant operator.

## Recommendations

It is recommended to always specialize on the properties of the operator as much as possible.

## Standard ODE Integrators

The standard ODE integrators will work on Non-autonomous linear ODE problems via an
automatic transformation to a first-order ODE. See the [ODE solvers](@ref ode_solve)
page for more details.

## Specialized OrdinaryDiffEq.jl Integrators

Unless otherwise specified, the OrdinaryDiffEq algorithms all come with a
3rd order Hermite polynomial interpolation. The algorithms denoted as having a
"free" interpolation means that no extra steps are required for the
interpolation. For the non-free higher order interpolating functions, the extra
steps are computed lazily (i.e. not during the solve).

Note that all of these methods are fixed timestep unless otherwise specified.

### Exponential Methods for Linear and Affine Problems

These methods require that ``A`` is constant.

- `LinearExponential` - Exact solution formula for linear, time-independent problems.

Options:

- `krylov` - symbol. One of
  - :off (default) - cache the operator beforehand. Requires `Matrix(A)` method
    defined for the operator `A`.
  - :simple - uses simple Krylov approximations with fixed subspace size `m`.
  - :adaptive - uses adaptive Krylov approximations with internal timestepping.
- `m` - integer, default: `30`. Controls the size of Krylov subsapce if
  `krylov=:simple`, and the initial subspace size if `krylov=:adaptive`.
- `iop` - integer, default: `0`. If not zero, determines the length of the incomplete
  orthogonalization procedure (IOP) [^1]. Note that if the linear operator/jacobian is hermitian,
  then the Lanczos algorithm will always be used and the IOP setting is ignored.
  
```julia
_A = [2 -1;-3 -5]/5
A = DiffEqArrayOperator(_A)
prob = ODEProblem(A, [1.0,-1.0], (1.0, 6.0))
sol = solve(prob, LinearExponential())
```


### State-Independent Solvers

These methods require ``A`` is only dependent on the independent variable, i.e. ``A(t)``.

- `MagnusMidpoint` - Second order Magnus Midpoint method.
- `MagnusLeapfrog`- Second order Magnus Leapfrog method.
- `MagnusGauss4` - Fourth order Magnus method approximated using a two stage Gauss quadrature.
- `MagnusGL4`- Fourth order Magnus method approximated using Gauss-Legendre quadrature.
- `MagnusNC6`- Sixth order Magnus method approximated using Newton-Cotes quadrature.
- `MagnusGL6`- Sixth order Magnus method approximated using Gauss-Legendre quadrature.
- `MagnusNC8`- Eighth order Magnus method approximated using Newton-Cotes quadrature.
- `MagnusGL8`- Eighth order Magnus method approximated using Gauss-Legendre quadrature.

Example:

```julia
function update_func(A,u,p,t)
    A[1,1] = cos(t)
    A[2,1] = sin(t)
    A[1,2] = -sin(t)
    A[2,2] = cos(t)
end
A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (1.0, 6.0))
sol = solve(prob,MagnusGL6(),dt=1/10)
```


### State-Dependent Solvers

These methods can be used when ``A`` is dependent on the state variables, i.e. ``A(u)``.


- `CayleyEuler` - First order method using Cayley transformations.
- `LieEuler` - First order Lie Euler method.
- `RKMK2` - Second order Runge–Kutta–Munthe-Kaas method.
- `RKMK4` - Fourth order Runge–Kutta–Munthe-Kaas method.
- `LieRK4` - Fourth order Lie Runge-Kutta method.
- `CG2` - Second order Crouch–Grossman method.
- `MagnusAdapt4` - Fourth Order Adaptive Magnus method.

Example:

```julia
function update_func(A,u,p,t)
    A[1,1] = 0
    A[2,1] = sin(u[1])
    A[1,2] = -1
    A[2,2] = 0
end
A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (0, 30.))
sol = solve(prob,LieRK4(),dt=1/4)
```

The above example solves a non-stiff Non-Autonomous Linear ODE
with a state dependent operator, using the `LieRK4` method.
Similarly, a stiff Non-Autonomous Linear ODE with state dependent
operators can be solved using specialized adaptive algorithms, like `MagnusAdapt4`. 

Example:

```julia
function update_func(A,u,p,t)
    A[1,1] = 0
    A[2,1] = 1
    A[1,2] = -2*(1 - cos(u[2]) - u[2]*sin(u[2]))
    A[2,2] = 0
end
A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (30, 150.))
sol = solve(prob,MagnusAdapt4())
```

# Time and State-Dependent Operators

These methods can be used when ``A`` is dependent on both time and state variables, i.e. ``A(u,t)``

- `CG3` - Third order Crouch-Grossman method.

[^1]: A description of IOP can be found in this [paper](https://doi.org/10.1016/j.jcp.2018.06.026).
