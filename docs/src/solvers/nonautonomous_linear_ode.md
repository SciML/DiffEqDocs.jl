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

Note that all of these methods are fixed timestep only.

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
