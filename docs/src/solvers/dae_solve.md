# DAE Solvers

## Recomended Methods

For medium to low accuracy DAEs in mass matrix form, the `Rodas4` and `Rodas42`
methods are good choices which will get good efficiency. The OrdinaryDiffEq.jl
methods are also the only methods which allow for Julia-defined number types.
For high accuracy (error `<1e-7`) on problems of `Vector{Float64}` defined in
mass matrix form, `radau` is an efficient method.

If the problem cannot be defined in mass matrix form, the recommended method for
performance is `IDA` from the Sundials.jl package if you are solving problems with
`Float64`. It's a very well-optimized method, and allows you to have a little bit of
control over the linear solver to better tailor it to your problem. A similar
algorithm is `daskr`. Which one is more efficient is problem-dependent.

# Full List of Methods

## OrdinaryDiffEq.jl

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

### Rosenbrock Methods

- `ROS3P` - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high
  accuracy on discretizations of nonlinear parabolic PDEs.
- `Rodas3` - 3rd order A-stable and stiffly stable Rosenbrock method.
- `Rodas4` - A 4th order A-stable stiffly stable Rosenbrock method with a
  stiff-aware 3rd order interpolant
- `Rodas42` - A 4th order A-stable stiffly stable Rosenbrock method with a
  stiff-aware 3rd order interpolant
- `Rodas4P` - A 4th order A-stable stiffly stable Rosenbrock method with a
  stiff-aware 3rd order interpolant. 4th order on linear parabolic problems
  and 3rd order accurate on nonlinear parabolic problems (as opposed to lower
  if not corrected).
- `Rodas5` - A 5th order A-stable stiffly stable Rosenbrock method. Currently has
  a Hermite interpolant because its stiff-aware 3rd order interpolant is not
  yet implemented.

### SDIRK Methods

SDIRK Methods

- `ImplicitEuler` - Stage order 1. A-B-L-stable. Adaptive
  timestepping through a divided differences estimate via memory. Strong-stability
  presurving (SSP).
- `ImplicitMidpoint` - Stage order 1. Symplectic. Good for when symplectic
  integration is required.

## Sundials.jl

- `IDA` - This is the IDA method from the Sundials.jl package.

Note that the constructors for the Sundials algorithms take a main argument:

- `linearsolver` - This is the linear solver which is used in the Newton iterations.
  The choices are:

  - `:Dense` - A dense linear solver
  - `:Band` - A solver specialized for banded Jacobians. If used, you must set the
    position of the upper and lower non-zero diagonals via `jac_upper` and
    `jac_lower`.
  - `BCG` - A Biconjugate gradient method.
  - `TFQMR` - A TFQMR method.

Example:

```julia
IDA() # Newton + Dense solver
IDA(linear_solver=:Band,jac_upper=3,jac_lower=3) # Banded solver with nonzero diagonals 3 up and 3 down
IDA(linear_solver=:BCG) # Biconjugate gradient method                                   
```

All of the additional options are available. The constructor is:

```julia
IDA(;linear_solver=:Dense,jac_upper=0,jac_lower=0,krylov_dim=0,
    max_order = 5,
    max_error_test_failures = 7,
    max_nonlinear_iters = 3,
    nonlinear_convergence_coefficient = 0.33,
    nonlinear_convergence_coefficient_ic = 0.0033,
    max_num_steps_ic = 5,
    max_num_jacs_ic = 4,
    max_num_iters_ic = 10,
    max_num_backs_ic = 100,
    use_linesearch_ic = true,
    max_convergence_failures = 10)
```

See [the Sundials manual](https://computation.llnl.gov/sites/default/files/public/ida_guide.pdf)
for details on the additional options.

## DASKR.jl

DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
Pkg.add("DASKR")
using DASKR
```

- `daskr` - This is a wrapper for the well-known DASKR algorithm.

All additional options are available. The constructor is:

```julia
function daskr(;linear_solver=:Dense,
                  jac_upper=0,jac_lower=0,max_order = 5,
                  non_negativity_enforcement = 0,
                  non_negativity_enforcement_array = nothing,
                  max_krylov_iters = nothing,
                  num_krylov_vectors = nothing,
                  max_number_krylov_restarts = 5,
                  krylov_convergence_test_constant = 0.05,
                  exclude_algebraic_errors = false)
```

Choices for the linear solver are:

- `:Dense`
- `:Banded`
- `:SPIGMR`, a Krylov method

## DASSL.jl

- `dassl` - A native Julia implementation of the DASSL algorithm.

## ODEInterfaceDiffEq.jl

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

- `seulex` - Extrapolation-algorithm based on the linear implicit Euler method.
- `radau` - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.
- `radau5` - Implicit Runge-Kutta method (Radau IIA) of order 5.
- `rodas` - Rosenbrock 4(3) method.
