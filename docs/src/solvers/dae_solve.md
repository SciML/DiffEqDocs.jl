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

## Full List of Methods

### OrdinaryDiffEq.jl

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

#### Rosenbrock Methods

- `ROS3P` - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high
  accuracy on discretizations of nonlinear parabolic PDEs.
- `Rodas3` - 3rd order A-stable and stiffly stable Rosenbrock method.
- `RosShamp4`- An A-stable 4th order Rosenbrock method.
- `Veldd4` - A 4th order D-stable Rosenbrock method.
- `Velds4` - A 4th order A-stable Rosenbrock method.
- `GRK4T` - An efficient 4th order Rosenbrock method.
- `GRK4A` - An A-stable 4th order Rosenbrock method. Essentially "anti-L-stable"
  but efficient.
- `Ros4LStab` - A 4th order L-stable Rosenbrock method.
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

#### Rosenbrock-W Methods

- `Rosenbrock23` - An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.
- `Rosenbrock32` - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildy stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.
- `RosenbrockW6S4OS` - A 4th order L-stable Rosenbrock-W method (fixed step only).
- `ROS34PW1a` - A 4th order L-stable Rosenbrock-W method.
- `ROS34PW1b` - A 4th order L-stable Rosenbrock-W method.
- `ROS34PW2` - A 4th order stiffy accurate Rosenbrock-W method for PDAEs.
- `ROS34PW3` - A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method.

#### SDIRK Methods

SDIRK Methods

- `ImplicitEuler` - Stage order 1. A-B-L-stable. Adaptive
  timestepping through a divided differences estimate via memory. Strong-stability
  presurving (SSP).
- `ImplicitMidpoint` - Stage order 1. Symplectic. Good for when symplectic
  integration is required.

### Sundials.jl

Note that this setup is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use Sundials.jl:

```julia
]add Sundials
using Sundials
```

- `IDA` - This is the IDA method from the Sundials.jl package.

Note that the constructors for the Sundials algorithms take a main argument:

- `linearsolver` - This is the linear solver which is used in the Newton iterations.
  The choices are:

  - `:Dense` - A dense linear solver.
  - `:Band` - A solver specialized for banded Jacobians. If used, you must set the
    position of the upper and lower non-zero diagonals via `jac_upper` and
    `jac_lower`.
  - `:GMRES` - A GMRES method. Recommended first choice Krylov method
  - `:BCG` - A Biconjugate gradient method.
  - `:PCG` - A preconditioned conjugate gradient method. Only for symmetric linear systems.
  - `:TFQMR` - A TFQMR method.
  - `:KLU` - A sparse factorization method. Requires that the user specifies a
    Jacobian. The Jacobian must be set as a sparse matrix in the `ODEProblem`
    type.

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
    max_convergence_failures = 10,
    init_all = false,
    prec = nothing, psetup = nothing, prec_side = 0)
```

See [the Sundials manual](https://computation.llnl.gov/sites/default/files/public/ida_guide.pdf)
for details on the additional options. The option `init_all` controls the initial condition
consistency routine. If the initial conditions are inconsistant (i.e. they do not satisfy the
implicit equation), `init_all=false` means that the algebraic variables and derivatives will
be modified in order to satisfy the DAE. If `init_all=true`, all initial conditions will be
modified to satify the DAE.

Note that here `prec` is a preconditioner function
`prec(z,r,p,t,y,fy,gamma,delta,lr)` where:

- `z`: the computed output vector
- `r`: the right-hand side vector of the linear system
- `p`: the parameters
- `t`: the current independent variable
- `du`: the current value of `f(u,p,t)`
- `gamma`: the `gamma` of `W = M - gamma*J`
- `delta`: the iterative method tolerance
- `lr`: a flag for whether `lr=1` (left) or `lr=2` (right)
  preconditioning

and `psetup` is the preconditioner setup function for pre-computing Jacobian
information. Where:

- `p`: the parameters
- `t`: the current independent variable
- `resid`: the current residual
- `u`: the current state
- `du`: the current `f(u,p,t)`
- `gamma`: the `gamma` of `W = M - gamma*J`

`psetup` is optional when `prec` is set.

### DASKR.jl

DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
]add DASKR
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

### DASSL.jl

- `dassl` - A native Julia implementation of the DASSL algorithm.

### ODEInterfaceDiffEq.jl

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

- `seulex` - Extrapolation-algorithm based on the linear implicit Euler method.
- `radau` - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.
- `radau5` - Implicit Runge-Kutta method (Radau IIA) of order 5.
- `rodas` - Rosenbrock 4(3) method.
