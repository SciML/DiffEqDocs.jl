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
- `KenCarp3`: An A-L stable stiffly-accurate 3rd order ESDIRK method
- `KenCarp4`: An A-L stable stiffly-accurate 4rd order ESDIRK method
- `KenCarp5`: An A-L stable stiffly-accurate 5rd order ESDIRK method

### Sundials.jl

- `ARKODE`: An additive Runge-Kutta method. Order between 3rd and 5th. For a list
  of available options, please see
  [its ODE solver page](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html#Sundials.jl-1)

## Semilinear ODE

The Semilinear ODE is a `SplitODEProblem` with one linear operator and one nonlinear function:

```math
\frac{du}{dt} =  Au + f(t,u)
```

See the documentation page for [DiffEqOperator](../../features/diffeq_operator.html) 
for details about how to define linear operators from a matrix or finite difference 
discretization of derivative operators.

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

- `GenericIIF1` - First order Implicit Integrating Factor method. Fixed timestepping only.
- `GenericIIF2` - Second order Implicit Integrating Factor method. Fixed timestepping only.
- `ETD1` - First order Exponential Time Differencing method. Not yet implemented.
- `ETD2` - Second order Exponential Time Differencing method. Not yet implemented.
- `LawsonEuler` - First order exponential Euler scheme. Fixed timestepping only.
- `NorsettEuler` - First order exponential-RK scheme. Fixed timestepping only.
- `ETDRK4` - 4th order exponential-RK scheme. Fixed timestepping only.

Note that the generic algorithms allow for a choice of `nlsolve`.

The methods need to compute the exponential of `A`, which could be expensive. There are 
two ways to speed up the integrator:

- For small systems that can fit `expm(dt*A)` in memory, use the in-place style, which 
  enables caching of the exponential operators to save time.

- For large systems, use Krylov-based versions of the methods which allow for lazy 
  calculation of `expm(dt*A)*v` and similar entities. To tell a solver to use Krylov 
  methods, pass `krylov=true` to its constructor. You can also manually set the size of the 
  Krylov subspace by setting the `m` parameter, which defaults to 30. For example
  
  ```julia
  LawsonEuler(krylob=true, m=50)
  ```
  
  constructs a Lawson-Euler method which uses a size-50 Krylov subspace. Note that `m` 
  only sets an upper bound to the Krylov subspace size. If a convergence criterion is met 
  (determined by the `reltol` of the integrator), "happy breakdown" will occur and the 
  Krylov subspace will only be constructed partially.
  
  Currently only the `LawsonEuler` and `NorsettEuler` methods support Krylov methods.
