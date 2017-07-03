# Refined ODE Solvers

`solve(prob::AbstractODEProblem,alg;kwargs)`

Solves the Refined ODE problems defined by `prob` using the algorithm `alg`.
If no algorithm is given, a default algorithm will be chosen.

## Special Forms

Many of the integrators in this category require special forms. For example,
sometimes an integrator may require that a certain argument is missing. Instead
of changing the function signature, keep the function signature but make sure
the function ignores the appropriate argument.

For example, one type of special form is the dynamical ODE:

```math
\frac{du}{dt} = f_1(t,v) \\
\frac{dv}{dt} = f_2(t,u) \\
```

This is a Partitioned ODE partitioned into two groups, so the functions should be
specified as `f1(t,x,v,dx)` and `f2(t,x,v,dx)` (in the inplace form). However,
this specification states that `f1` would be independent of `x`, and `f2` should
be independent of `v`. Following the requirements for the integrator is required
to achieve the suggested accuracy.

## Note About OrdinaryDiffEq.jl

Unless otherwise specified, the OrdinaryDiffEq algorithms all come with a
3rd order Hermite polynomial interpolation. The algorithms denoted as having a
"free" interpolation means that no extra steps are required for the
interpolation. For the non-free higher order interpolating functions, the extra
steps are computed lazily (i.e. not during the solve).

# Functional Forms

## Dynamical ODE

These algorithms require a Partitioned ODE of the form:

```math
\frac{du}{dt} = f_1(v) \\
\frac{dv}{dt} = f_2(t,u) \\
```
This is a Partitioned ODE partitioned into two groups, so the functions should be
specified as `f1(t,u,v,dx)` and `f2(t,u,v,dv)` (in the inplace form), where `f1`
is independent of `t` and `u`, and `f2` is independent of `v`. This includes
discretizations arising from `SecondOrderODEProblem`s where the velocity is not
used in the acceleration function, and Hamiltonians where the potential is (or
can be) time-dependent but the kinetic energy is only dependent on `v`.

Note that some methods assume that the integral of `f1` is a quadratic form. That
means that `f1=v'*M*v`, i.e. ``\int f1 = 1/2 m v^2``, giving `du = v`. This is
equivalent to saying that the kinetic energy is related to ``v^2``. The methods
which require this assumption will lose accuracy if this assumption is violated.
Methods listed below make note of this requirement with "Requires quadratic
kinetic energy".

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

- `SymplecticEuler`: First order explicit symplectic integrator
- `VelocityVerlet`: 2nd order explicit symplectic integrator.
- `VerletLeapfrog`: 2nd order explicit symplectic integrator.
- `PseudoVerletLeapfrog`: 2nd order explicit symplectic integrator.
- `McAte2`: Optimized efficiency 2nd order explicit symplectic integrator.
- `Ruth3`: 3rd order explicit symplectic integrator.
- `McAte3`: Optimized efficiency 3rd order explicit symplectic integrator.
- `CandyRoz4`: 4th order explicit symplectic integrator.
- `McAte4`: 4th order explicit symplectic integrator. Requires quadratic
  kinetic energy.
- `CalvoSanz4`: Optimized efficiency 4th order explicit symplectic integrator.
- `McAte42`: 4th order explicit symplectic integrator.
- `McAte5`: Optimized efficiency 5th order explicit symplectic integrator.
  Requires quadratic kinetic energy
- `Yoshida6`: 6th order explicit symplectic integrator.
- `KahanLi6`: Optimized efficiency 6th order explicit symplectic integrator.
- `McAte8`: 8th order explicit symplectic integrator.
- `KahanLi8`: Optimized efficiency 8th order explicit symplectic integrator.
- `SofSpa10`: 10th order explicit symplectic integrator.

### Recommendations

Higher order algorithms are the most efficient when higher accuracy is needed,
and when less accuracy is needed lower order methods do better. Optimized efficiency
methods take more steps and thus have more force calculations for the same order,
but have smaller error. Thus the "optimized efficiency" algorithms are recommended
if your force calculation is not too sufficiency large, while the other methods are
recommend when force calculations are really large (for example, like in MD simulations
`VelocityVerlet` is very popular since it only requires one force calculation
per timestep).

## Implicit-Explicit (IMEX) ODE

The Implicit-Explicit (IMEX) ODE is a split `ODEProblem` with two functions:

```math
\frac{du}{dt} =  f_1(t,u) + f_2(t,u)
```

where the first function is the stiff part and the second function is the non-stiff
part (implicit integration on `f1`, explicit integration on `f2`).

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

- `SplitEuler`: 1st order fully explicit method. Used for testing accuracy
  of splits.

### Sundials.jl

- `ARKODE`: An additive Runge-Kutta method. Not yet implemented.

## Linear-Nonlinear (LNL) ODE

The Linear-Nonlinear (LNL) ODE is a split `ODEProblem` with two functions:

```math
\frac{du}{dt} =  f_1(t,u) + f_2(t,u)
```

where the first function is a `DiffEqOperator` and the second function is the
non-stiff part (implicit integration on `f1`, explicit integration on `f2`).

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

- `IIF1` - First order Implicit Integrating Factor method. Not yet implemented.
- `IIF2` - Second order Implicit Integrating Factor method. Not yet implemented.
- `ETD1` - First order Exponential Time Differencing method. Not yet implemented.
- `ETD2` - Second order Exponential Time Differencing method. Not yet implemented.
- `ExpEuler` - First order exponential Euler scheme. Not yet implemented.
- `NorsettEuler` - First order exponential-RK scheme. Not yet implemented.
