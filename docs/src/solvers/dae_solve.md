# DAE Solvers

## Recomended Methods

The recommended method for performance is `IDA` from the Sundials.jl package if you
are solving problems with `Float64`. It's a very well-optimized method, and allows
you to have a little bit of control over the linear solver to better tailor it
to your problem. A similar algorithm is `daskr`. Which one is more efficient
is problem-dependent.

If your problem requires special Julia types like arbitrary precision numbers,
then `dassl` is the method of choice.

## Special Keyword Arguments

## Implemented Solvers

### Sundials.jl

- `IDA` - This is the IDA method from the Sundials.jl package.

Note that the constructors for the Sundials algorithms take a main argument:

- `linearsolver` - This is the linear solver which is used in the Newton iterations.
  The choices are:

  - `:Dense` - A dense linear solver
  - `:Band` - A solver specialized for banded Jacobians. If used, you must set the
    position of the upper and lower non-zero diagonals via `jac_upper` and
    `jac_lower`.
  - `:Diagonal` - This method is specialized for diagonal Jacobians.
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

### DASKR.jl

DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
Pkg.add("DASKR")
using DASKR
```

- `daskr` - This is a wrapper for the well-known DASKR algorithm.

### DASSL.jl

- `dassl` - A native Julia implementation of the DASSL algorithm.
