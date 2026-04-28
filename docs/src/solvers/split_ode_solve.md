# [Split ODE Solvers](@id split_ode_solve)

The solvers which are available for a `SplitODEProblem` depend on the input
linearity and number of components. Each solver has functional form
(or many) that it allows.

## Packages

The solvers on this page are distributed across the packages below. Add the package(s) you need to your environment.

| Package | Description |
|---|---|
| `OrdinaryDiffEqLowOrderRK` | Low-order explicit Runge-Kutta methods (BS3, DP5, RK4, Heun, Euler, OwrenZen, etc.). |
| `OrdinaryDiffEqSDIRK` | SDIRK / ESDIRK methods for stiff problems (KenCarp3/4/47/58, TRBDF2, ImplicitEuler, Kvaerno). |
| `OrdinaryDiffEqBDF` | BDF / NDF multistep (FBDF, QNDF, ABDF2, SBDF) and implicit-DAE forms (DFBDF, DImplicitEuler). |
| `OrdinaryDiffEqIMEXMultistep` | IMEX multistep methods for split ODEs (CNAB2, CNLF2). |
| `OrdinaryDiffEqExponentialRK` | Exponential RK and exponential propagation methods (LawsonEuler, ETDRK4, EPIRK, Exprb). |
| `Sundials` | Wrappers for the SUNDIALS C library: `CVODE_BDF`, `CVODE_Adams`, `IDA`, `ARKODE`. |


## Implicit-Explicit (IMEX) ODE

The Implicit-Explicit (IMEX) ODE is a `SplitODEProblem` with two functions:

```math
\frac{du}{dt} =  f_1(t,u) + f_2(t,u)
```

where the first function is the stiff part and the second function is the non-stiff
part (implicit integration on `f1`, explicit integration on `f2`).

### Recommended Methods

The recommended method in most cases is `KenCarp4`. In cases of extreme stiffness
or for high tolerances, `KenCarp3` can be a good choice. The `ARKODE` methods
are generally inefficient and diverge unless the options are tweaked to match
the problem, though for large enough PDEs the `ARKODE` method with
`linear_solver=:GMRES` is a good choice.

### OrdinaryDiffEq.jl

!!! note "v8: import from the IMEX OrdinaryDiffEq sublibs"

    None of the IMEX solvers below are in OrdinaryDiffEq's default re-export
    set under v7.

    ```julia
    using OrdinaryDiffEqLowOrderRK    # SplitEuler
    using OrdinaryDiffEqBDF           # IMEXEuler, IMEXEulerARK, SBDF2, SBDF3, SBDF4, MEBDF2
    using OrdinaryDiffEqIMEXMultistep # CNAB2, CNLF2
    using OrdinaryDiffEqSDIRK         # KenCarp3, KenCarp4, KenCarp47, KenCarp5, KenCarp58
    ```

  - `OrdinaryDiffEqLowOrderRK.SplitEuler`: 1st order fully explicit method. Used for testing accuracy
    of splits.
  - `OrdinaryDiffEqBDF.IMEXEuler` : 1st order explicit Euler mixed with implicit Euler. Fixed time
    step only.
  - `OrdinaryDiffEqIMEXMultistep.CNAB2`: Crank-Nicolson Adams Bashforth Order 2. Fixed time step only.
  - `OrdinaryDiffEqIMEXMultistep.CNLF2`: Crank-Nicolson Leapfrog of Order 2. Fixed time step only.
  - `OrdinaryDiffEqBDF.SBDF2` : 2nd order IMEX BDF method. Fixed time step only.
  - `OrdinaryDiffEqBDF.SBDF3` : 3rd order IMEX BDF method. Fixed time step only. In development.
  - `OrdinaryDiffEqBDF.SBDF4` : 4th order IMEX BDF method. Fixed time step only. In development.
  - `OrdinaryDiffEqSDIRK.KenCarp3`: An A-L stable stiffly-accurate 3rd order ESDIRK method.
  - `OrdinaryDiffEqSDIRK.KenCarp4`: An A-L stable stiffly-accurate 4th order ESDIRK method.
  - `OrdinaryDiffEqSDIRK.KenCarp47` - An A-L stable stiffly-accurate 4th order seven-stage ESDIRK method with splitting
  - `OrdinaryDiffEqSDIRK.KenCarp5`: An A-L stable stiffly-accurate 5th order ESDIRK method.
  - `OrdinaryDiffEqSDIRK.KenCarp58` - An A-L stable stiffly-accurate 5th order eight-stage ESDIRK method with splitting

### Sundials.jl

  - `Sundials.ARKODE`: An additive Runge-Kutta method. Order between 3rd and 5th. For a list
    of available options, please see
    [its ODE solver page](https://docs.sciml.ai/DiffEqDocs/dev/api/sundials/).

## Semilinear ODE

The Semilinear ODE is a `SplitODEProblem` with one linear operator and one nonlinear function:

```math
\frac{du}{dt} =  Au + f(t,u)
```

See the documentation page for [SciMLOperators](https://docs.sciml.ai/SciMLOperators/stable/)
for details about how to define linear operators from a matrix or finite difference
discretization of derivative operators.

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

!!! note "v8: import from `OrdinaryDiffEqExponentialRK`"

    The exponential RK / EPIRK families below all live in
    `OrdinaryDiffEqExponentialRK` (a sublib of OrdinaryDiffEq). Bring them in
    with `using OrdinaryDiffEqExponentialRK`.

  - `OrdinaryDiffEqExponentialRK.LawsonEuler` - First order exponential Euler scheme. Fixed timestepping only.
  - `OrdinaryDiffEqExponentialRK.NorsettEuler` - First order exponential-RK scheme. Fixed timestepping only. Alias: `ETD1`.
  - `OrdinaryDiffEqExponentialRK.ETD2` - Second order Exponential Time Differencing method (in development). Fixed timestepping only. Doesn't support Krylov approximation.
  - `OrdinaryDiffEqExponentialRK.ETDRK2` - 2nd order exponential-RK scheme. Fixed timestepping only.
  - `OrdinaryDiffEqExponentialRK.ETDRK3` - 3rd order exponential-RK scheme. Fixed timestepping only.
  - `OrdinaryDiffEqExponentialRK.ETDRK4` - 4th order exponential-RK scheme. Fixed timestepping only.
  - `OrdinaryDiffEqExponentialRK.HochOst4` - 4th order exponential-RK scheme with stiff order 4. Fixed
    timestepping only.

Note that the generic algorithms `GenericIIF1` and `GenericIIF2` allow for a choice of `nlsolve`.

By default, the exponential methods cache matrix functions such as `exp(dt*A)` to accelerate
the time stepping for small systems. For large systems, using Krylov-based versions of the
methods can allow for lazy calculation of `exp(dt*A)*v` and similar entities, and thus improve
performance.

To tell a solver to use Krylov methods, pass `krylov=true` to its constructor. You
can also manually set the size of the Krylov subspace by setting the `m` parameter, which
defaults to 30. For example

```julia
LawsonEuler(krylov = true, m = 50)
```

constructs a Lawson-Euler method, which uses a size-50 Krylov subspace. Note that `m`
only sets an upper bound to the Krylov subspace size. If a convergence criterion is met
(determined by the `reltol` of the integrator), “happy breakdown” will occur and the
Krylov subspace will only be constructed partially.

For more advanced control over the Krylov algorithms, you can change the length of the
incomplete orthogonalization procedure (IOP) [^1] by setting the `iop` parameter in the
constructor. By default, IOP is turned off and full Arnoldi iteration is used. Note that
if the linear operator is hermitian, then the Lanczos algorithm will always be used and
IOP setting is ignored.

[^1]: Koskela, A. (2015). Approximating the matrix exponential of an advection-diffusion operator using the incomplete orthogonalization method. In Numerical Mathematics and Advanced Applications-ENUMATH 2013 (pp. 345-353). Springer, Cham.
