# Dynamical, Hamiltonian, and 2nd Order ODE Solvers

Dynamical ODEs, such as those arising from Hamiltonians or second
order ordinary differential equations, give rise to a special structure
that can be specialized on in the solver for more efficiency.
These algorithms require an ODE defined in the following ways:

```julia
DynamicalODEProblem{isinplace}(f1, f2, v0, u0, tspan, p = NullParameters(); kwargs...)
SecondOrderODEProblem{isinplace}(f, du0, u0, tspan, p = NullParameters(); kwargs...)
HamiltonianProblem{T}(H, p0, q0, tspan, p = NullParameters(); kwargs...)
```

These correspond to partitioned equations of motion:

```math
\frac{dv}{dt} = f_1(t,u) \\
\frac{du}{dt} = f_2(v) \\
```

or, for `SecondOrderODEProblem`,

```math
\frac{d^2u}{dt^2} = f(u,p,t)
```

The functions should be specified as `f1(dv,v,u,p,t)` and `f2(du,v,u,p,t)`
(in the inplace form), where `f1` is independent of `v` (unless
specified by the solver), and `f2` is independent of `t` and `u`. This includes
discretizations arising from `SecondOrderODEProblem`s where the velocity is not
used in the acceleration function, and Hamiltonians where the potential is
(or can be) time-dependent, but the kinetic energy is only dependent on `v`.

Note that some methods assume that the integral of `f2` is a quadratic form. That
means that `f2=v'*M*v`, i.e. ``\int f_2 = \frac{1}{2} m v^2``, giving `du = v`. This is
equivalent to saying that the kinetic energy is related to ``v^2``. The methods
which require this assumption will lose accuracy if this assumption is violated.
Methods listed below make note of this requirement with "Requires quadratic
kinetic energy".

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package | Description |
|---|---|
| `OrdinaryDiffEqRKN` | Runge-Kutta-Nystrom methods (DPRKN6/8/12, ERKN4/5/7, IRKN3/4) for second-order ODEs. |
| `OrdinaryDiffEqSymplecticRK` | Symplectic integrators (KahanLi6/8, McAte, VelocityVerlet, Yoshida6) for Hamiltonian systems. |
| `GeometricIntegratorsDiffEq` | Wrappers for GeometricIntegrators.jl (Gauss, Lobatto, Radau, Symplectic methods). |


## Recommendations

When energy conservation is required, use a symplectic method. Otherwise, the
Runge-Kutta-Nyström methods will be more efficient. Energy is mostly conserved
by Runge-Kutta-Nyström methods, but is not conserved for long-time integrations.
Thus, it is suggested that for shorter integrations you use Runge-Kutta-Nyström
methods as well.

As a go-to method for efficiency, `DPRKN6` is a good choice. `DPRKN12` is a good
choice when high accuracy, like `tol<1e-10` is necessary. However, `DPRKN6` is
the only Runge-Kutta-Nyström method with a higher order interpolant (all default
to order 3 Hermite, whereas `DPRKN6` is order 6th interpolant) and thus in cases
where interpolation matters (ex: event handling) one should use `DPRKN6`. For
very smooth problems with expensive acceleration function evaluations, `IRKN4`
can be a good choice as it minimizes the number of evaluations.

For symplectic methods, higher order algorithms are the most efficient when higher
accuracy is needed, and when less accuracy is needed lower order methods do better.
Optimized efficiency methods take more steps and thus have more force calculations
for the same order, but have smaller error. Thus, the “optimized efficiency”
algorithms are recommended if your force calculation is not too sufficiency large,
while the other methods are recommended when force calculations are really large
(for example, like in MD simulations `VelocityVerlet` is very popular since it only
requires one force calculation per timestep). A good go-to method would be `McAte5`,
and a good high order choice is `KahanLi8`.

## Standard ODE Integrators

The standard ODE integrators will work on Dynamical ODE problems via an automatic
transformation to a first-order ODE. See the [ODE solvers](@ref ode_solve)
page for more details.

## Specialized OrdinaryDiffEq.jl Integrators

Unless otherwise specified, the OrdinaryDiffEq algorithms all come with a
3rd order Hermite polynomial interpolation. The algorithms denoted as having a
“free” interpolation means that no extra steps are required for the
interpolation. For the non-free higher order interpolating functions, the extra
steps are computed lazily (i.e. not during the solve).

!!! note "v8: import from the dynamical OrdinaryDiffEq sublibs"

    None of the specialized solvers below are in OrdinaryDiffEq's default
    re-export set under v7. Bring them in directly:

    ```julia
    using OrdinaryDiffEqRKN          # Nystrom*, IRKN*, ERKN*, DPRKN*
    using OrdinaryDiffEqSymplecticRK # SymplecticEuler, VelocityVerlet, ...
                                     # McAte*, KahanLi6, KahanLi8, SofSpa10
    ```

### Runge-Kutta-Nyström Integrators

  - `OrdinaryDiffEqRKN.Nystrom4`: 4th order explicit Runge-Kutta-Nyström method. Allows acceleration
    to depend on velocity. Fixed timestep only.
  - `OrdinaryDiffEqRKN.IRKN3`: 4th order explicit two-step Runge-Kutta-Nyström method. Fixed
    timestep only.
  - `OrdinaryDiffEqRKN.IRKN4`: 4th order explicit two-step Runge-Kutta-Nyström method. Can be more
    efficient for smooth problems. Fixed timestep only.
  - `OrdinaryDiffEqRKN.ERKN4`: 4th order Runge-Kutta-Nyström method which integrates the periodic
    properties of the harmonic oscillator exactly. Gets extra efficiency on periodic
    problems.
  - `OrdinaryDiffEqRKN.ERKN5`: 5th order Runge-Kutta-Nyström method which integrates the periodic
    properties of the harmonic oscillator exactly. Gets extra efficiency on periodic
    problems.
  - `OrdinaryDiffEqRKN.ERKN7`: 7th order Runge-Kutta-Nyström method which integrates the periodic
    properties of the harmonic oscillator exactly. Gets extra efficiency on periodic
    problems.
  - `OrdinaryDiffEqRKN.Nystrom4VelocityIndependent`: 4th order explicit Runge-Kutta-Nyström method.
    Fixed timestep only.
  - `OrdinaryDiffEqRKN.Nystrom5VelocityIndependent`: 5th order explicit Runge-Kutta-Nyström method.
    Fixed timestep only.
  - `OrdinaryDiffEqRKN.DPRKN4`: 4th order explicit adaptive Runge-Kutta-Nyström method.
  - `OrdinaryDiffEqRKN.DPRKN5`: 5th order explicit adaptive Runge-Kutta-Nyström method.
  - `OrdinaryDiffEqRKN.DPRKN6`: 6th order explicit adaptive Runge-Kutta-Nyström method. Free 6th
    order interpolant.
  - `OrdinaryDiffEqRKN.DPRKN6FM`: 6th order explicit adaptive Runge-Kutta-Nyström method.
  - `OrdinaryDiffEqRKN.DPRKN8`: 8th order explicit adaptive Runge-Kutta-Nyström method.
  - `OrdinaryDiffEqRKN.DPRKN12`: 12th order explicit adaptive Runge-Kutta-Nyström method.

### Symplectic Integrators

Note that all symplectic integrators are fixed timestep only.

  - `OrdinaryDiffEqSymplecticRK.SymplecticEuler`: First order explicit symplectic integrator
  - `OrdinaryDiffEqSymplecticRK.VelocityVerlet`: 2nd order explicit symplectic integrator. Requires `f_2(t,u) = v`, i.e.
    a second order ODE.
  - `OrdinaryDiffEqSymplecticRK.VerletLeapfrog`: 2nd order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.PseudoVerletLeapfrog`: 2nd order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.McAte2`: Optimized efficiency 2nd order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.Ruth3`: 3rd order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.McAte3`: Optimized efficiency 3rd order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.CandyRoz4`: 4th order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.McAte4`: 4th order explicit symplectic integrator. Requires quadratic
    kinetic energy.
  - `OrdinaryDiffEqSymplecticRK.CalvoSanz4`: Optimized efficiency 4th order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.McAte42`: 4th order explicit symplectic integrator. (Broken)
  - `OrdinaryDiffEqSymplecticRK.McAte5`: Optimized efficiency 5th order explicit symplectic integrator.
    Requires quadratic kinetic energy.
  - `OrdinaryDiffEqSymplecticRK.Yoshida6`: 6th order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.KahanLi6`: Optimized efficiency 6th order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.McAte8`: 8th order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.KahanLi8`: Optimized efficiency 8th order explicit symplectic integrator.
  - `OrdinaryDiffEqSymplecticRK.SofSpa10`: 10th order explicit symplectic integrator.

### GeometricIntegrators.jl

GeometricIntegrators.jl is a set of fixed timestep algorithms written in Julia.
Note that this setup is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use
GeometricIntegratorsDiffEq.jl:

```julia
Pkg.clone("https://github.com/SciML/GeometricIntegratorsDiffEq.jl")
import GeometricIntegratorsDiffEq
```

  - `GeometricIntegratorsDiffEq.GISymplecticEulerA` - First order explicit symplectic Euler A
  - `GeometricIntegratorsDiffEq.GISymplecticEulerB` - First order explicit symplectic Euler B
  - `GeometricIntegratorsDiffEq.GILobattoIIIAIIIB(n)` - Nth order Gauss-Labatto-IIIA-IIIB
  - `GeometricIntegratorsDiffEq.GILobattoIIIBIIIA(n)` - Nth order Gauss-Labatto-IIIB-IIIA
