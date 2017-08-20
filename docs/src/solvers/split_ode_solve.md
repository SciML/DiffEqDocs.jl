# Split ODE Solvers

The solvers which are available for a `SplitODEProblem` depend on the input
linearity and number of components. Each solver has functional form
(or many) that it allows.

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

## Semilinear ODE

The Semilinear ODE is a split `ODEProblem` with two functions:

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
