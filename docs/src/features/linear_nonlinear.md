# Specifying (Non)Linear Solvers

One of the key features of DifferentialEquations.jl is its flexibility. Keeping
with this trend, many of the native Julia solvers provided by DifferentialEquations.jl
allow you to choose the method for linear and nonlinear solving. This section
details how to make that choice.

## Linear Solvers: `linsolve` Specification

For differential equation integrators which use linear solvers, an argument
to the method `linsolve` determines the linear solver which is used. The signature
is:

```julia
linsolve! = linsolve(Val{:init},f,x;kwargs...)
linsolve!(x,A,b,matrix_updated=false;kwargs...)
```

This is an in-place function which updates `x` by solving `Ax=b`. The user should
specify the function `linsolve(Val{:init},f,x)` which returns a `linsolve!` function.
The setting `matrix_updated` determines whether the matrix `A` has changed from the
last call. This can be used to smartly cache factorizations.

Note that `linsolve!` needs to accept splatted keyword arguments. The possible arguments
passed to the linear solver are as follows:

- `Pl`, a pre-specified left preconditioner which utilizes the internal adaptive norm estimates
- `Pr`, a pre-specified right preconditioner which utilizes the internal adaptive norm estimates
- `tol`, a linear solver tolerance specified from the ODE solver's implicit handling

### Pre-Built Linear Solver Choices

The following choices of pre-built linear solvers exist:

- DefaultLinSolve
- LinSolveFactorize
- LinSolveGPUFactorize
- LinSolveGMRES
- LinSolveCG
- LinSolveBiCGStabl
- LinSolveChebyshev
- LinSolveMINRES
- LinSolveIterativeSolvers

### DefaultLinSolve

The default linear solver is `DefaultLinSolve`. This method is adaptive, and
automatically chooses an LU factorization choose for dense and sparse
arrays, and is compatible with GPU-based arrays. When the Jacobian is an
`AbstractDiffEqOperator`, i.e. is matrix-free, `DefaultLinSolve` defaults to
using a `gmres` iterative solver.

### Basic linsolve method choice: Factorization by LinSolveFactorize

The easiest way to specify a `linsolve` is by a `factorization` function which
generates a type on which `\` (or `A_ldiv_B!`) is called.  This is done through
the helper function `LinSolveFactorize` which makes the appropriate function.
For example, the  `Rosenbrock23` takes in a `linsolve` function, which we can
choose to be a  QR-factorization by:

```julia
Rosenbrock23(linsolve=LinSolveFactorize(qrfact!))
```

LinSolveFactorize takes in a function which returns an object that can `\`.
Direct methods like `qrfact!` will automatically cache the factorization,
making it efficient for small dense problems.

However, for large sparse problems, you can let `\` be an iterative method. For
example, using PETSc.jl, we can define our factorization function to be:

```julia
linsolve = LinSolveFactorize((A) -> KSP(A, ksp_type="gmres", ksp_rtol=1e-6))
```

This function creates a `KSP` type which makes `\` perform the GMRES iterative
method provided by PETSc.jl. Thus if we pass this function into the algorithm
as the factorization method, all internal linear solves will happen by PETSc.jl.

### GPU offloading of factorization with LinSolveGPUFactorize

If one has a problem with a sufficiently large Jacobian (~100x100) and a
sufficiently powerful GPU, it can make sense to offload the factorization
and backpropogation steps to the GPU. For this, the `LinSolveGPUFactorize`
linear solver is provided. It works similarly to `LinSolveFactorize`, but
the matrix is automatically sent to the GPU as a `CuArray` and the `ldiv!`
is performed against a CUDA QR factorization of the matrix.

Note that this method requires that you have done `using CuArrays` in your
script. A working installation of CuArrays.jl is required, which requires
an installation of CUDA Toolkit.

### IterativeSolvers.jl-Based Methods

The signature for `LinSolveIterativeSolvers` is:

```julia
LinSolveIterativeSolvers(generate_iterator,args...;
                         Pl=IterativeSolvers.Identity(),
                         Pr=IterativeSolvers.Identity(),
                         kwargs...)
```

where `Pl` is the left preconditioner, `Pr` is the right preconditioner, and
the other `args...` and `kwargs...` are passed into the iterative solver
chosen in `generate_iterator` which designates the construction of an iterator
from IterativeSolvers.jl. For example, using `gmres_iterable!` would make a
version that uses `IterativeSolvers.gmres`. The following are aliases to common
choices:

- LinSolveGMRES -- GMRES
- LinSolveCG -- CG (Conjugate Gradient)
- LinSolveBiCGStabl -- BiCGStabl Stabilized Bi-Conjugate Gradient
- LinSolveChebyshev -- Chebyshev
- LinSolveMINRES -- MINRES

which all have the same arguments as `LinSolveIterativeSolvers` except with
`generate_iterator` pre-specified.

### Implementing Your Own LinSolve: How LinSolveFactorize Was Created

In order to make your own `linsolve` functions, let's look at how the `LinSolveFactorize`
function is created. For example, for an LU-Factorization, we would like to use
`lufact!` to do our linear solving. We can directly write this as:

```julia
function linsolve!(::Type{Val{:init}},f,u0,kwargs...)
  function _linsolve!(x,A,b,update_matrix=false,kwargs...)
    _A = lufact!(A)
    ldiv!(x,_A,b)
  end
end
```

This initialization function returns a linear solving function
that always computes the LU-factorization and then does the solving.
This method works fine and you can pass it to the methods like

```julia
Rosenbrock23(linsolve=linsolve!)
```

and it will work, but this method does not cache `_A`, the factorization. This
means that, even if `A` has not changed, it will re-factorize the matrix.

To change this, we can instead create a call-overloaded type. The generalized form
of this is:

```julia
mutable struct LinSolveFactorize{F}
  factorization::F
  A
end
LinSolveFactorize(factorization) = LinSolveFactorize(factorization,nothing)
function (p::LinSolveFactorize)(x,A,b,matrix_updated=false)
  if matrix_updated
    p.A = p.factorization(A)
  end
  A_ldiv_B!(x,p.A,b)
end
function (p::LinSolveFactorize)(::Type{Val{:init}},f,u0_prototype)
  LinSolveFactorize(p.factorization,nothing)
end
linsolve = LinSolveFactorize(lufact!)
```

`LinSolveFactorize` is a type which holds the factorization method and the pre-factorized
matrix. When `linsolve` is passed to the ODE/SDE/etc. solver, it will use the function
`linsolve(Val{:init},f,u0_prototype)` to create a `LinSolveFactorize` object which holds
the factorization method and a cache for holding a factorized matrix. Then

```julia
function (p::LinSolveFactorize)(x,A,b,matrix_updated=false)
  if matrix_updated
    p.A = p.factorization(A)
  end
  A_ldiv_B!(x,p.A,b)
end
```

is what's used in the solver's internal loop. If `matrix_updated` is true, it
will re-compute the factorization. Otherwise it just solves the linear system
with the cached factorization. This general idea of using a call-overloaded
type can be employed to do many other things.
