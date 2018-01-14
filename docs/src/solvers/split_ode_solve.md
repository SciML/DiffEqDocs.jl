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

The Semilinear ODE is a split `ODEProblem` with one linear operator and one function:

```math
\frac{du}{dt} =  Au + f(t,u)
```

where the first function is a constant (not time dependent)`AbstractDiffEqOperator`
and the second part is a (nonlinear) function.
[../../features/diffeq_operator.html](See the DiffEqOperator page for details).

The appropriate algorithms for this form are:

### OrdinaryDiffEq.jl

- `GenericIIF1` - First order Implicit Integrating Factor method. Fixed timestepping only.
- `GenericIIF2` - Second order Implicit Integrating Factor method. Fixed timestepping only.
- `ETD1` - First order Exponential Time Differencing method. Not yet implemented.
- `ETD2` - Second order Exponential Time Differencing method. Not yet implemented.
- `LawsonEuler` - First order exponential Euler scheme. Fixed timestepping only.
- `NorsettEuler` - First order exponential-RK scheme. Fixed timestepping only.

Note that the generic algorithms allow for a choice of `nlsolve`.
