# Split ODE Solvers

The solvers which are available for a `SplitODEProblem` depend on the input
linearity and number of components. Each solver has functional form
(or many) that it allows.

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

- `SplitEuler`: 1st order fully explicit method. Used for testing accuracy
  of splits.
- `IMEXEuler` : 1st order explicit Euler mixed with implicit Euler. Fixed time
  step only.
- `CNAB2`: Crank-Nicholson Adams Bashforth Order 2. Fixed time step only.
- `CNLF`: Crank-Nicholson Leapfrog of Order 2. Fixed time step only.
- `SBDF2` : 2nd order IMEX BDF method. Fixed time step only.
- `SBDF3` : 3rd order IMEX BDF method. Fixed time step only. In development.
- `SBDF4` : 4th order IMEX BDF method. Fixed time step only. In development.
- `KenCarp3`: An A-L stable stiffly-accurate 3rd order ESDIRK method.
- `KenCarp4`: An A-L stable stiffly-accurate 4rd order ESDIRK method.
- `KenCarp5`: An A-L stable stiffly-accurate 5rd order ESDIRK method.

### Sundials.jl

- `ARKODE`: An additive Runge-Kutta method. Order between 3rd and 5th. For a list
  of available options, please see
  [its ODE solver page](http://docs.juliadiffeq.org/latest/solvers/ode_solve#Sundials.jl-1)

## Semilinear ODE

The Semilinear ODE is a `SplitODEProblem` with one linear operator and one nonlinear function:

```math
\frac{du}{dt} =  Au + f(t,u)
```

See the documentation page for [DiffEqOperator](../../../features/diffeq_operator)
for details about how to define linear operators from a matrix or finite difference
discretization of derivative operators.

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

- `GenericIIF1` - First order Implicit Integrating Factor method. Fixed timestepping only. Doesn't support Krylov approximation.
- `GenericIIF2` - Second order Implicit Integrating Factor method. Fixed timestepping only. Doesn't support Krylov approximation.
- `LawsonEuler` - First order exponential Euler scheme. Fixed timestepping only.
- `NorsettEuler` - First order exponential-RK scheme. Fixed timestepping only. Alias: `ETD1`.
- `ETD2` - Second order Exponential Time Differencing method (in development). Fixed timestepping only. Doesn't support Krylov approximation.
- `ETDRK2` - 2nd order exponential-RK scheme. Fixed timestepping only.
- `ETDRK3` - 3rd order exponential-RK scheme. Fixed timestepping only.
- `ETDRK4` - 4th order exponential-RK scheme. Fixed timestepping only.
- `HochOst4` - 4th order exponential-RK scheme with stiff order 4. Fixed
  timestepping only.
- `Exprb32` - 3rd order adaptive Exponential Rosenbrock scheme (in development).
- `Exprb43` - 4th order adaptive Exponential Rosenbrock scheme (in development).

Note that the generic algorithms `GenericIIF1` and `GenericIIF2` allow for a choice of `nlsolve`.

By default, the exponential methods cache matrix functions such as `exp(dt*A)` to accelerate
the time stepping for small systems. For large systems, using Krylov-based versions of the
methods can allow for lazy calculation of `exp(dt*A)*v` and similar entities, and thus improve
performance.

To tell a solver to use Krylov methods, pass `krylov=true` to its constructor. You
can also manually set the size of the Krylov subspace by setting the `m` parameter, which
defaults to 30. For example

```julia
LawsonEuler(krylov=true, m=50)
```

constructs a Lawson-Euler method which uses a size-50 Krylov subspace. Note that `m`
only sets an upper bound to the Krylov subspace size. If a convergence criterion is met
(determined by the `reltol` of the integrator), "happy breakdown" will occur and the
Krylov subspace will only be constructed partially.

For more advanced control over the Krylov algorithms, you can change the length of the
incomplete orthogonalization procedure (IOP) [^1] by setting the `iop` parameter in the
constructor. By default, IOP is turned off and full Arnoldi iteration is used. Note that
if the linear operator is hermitian, then the Lanczos algorithm will always be used and
IOP setting is ignored.

[^1]: Koskela, A. (2015). Approximating the matrix exponential of an advection-diffusion operator using the incomplete orthogonalization method. In Numerical Mathematics and Advanced Applications-ENUMATH 2013 (pp. 345-353). Springer, Cham.
