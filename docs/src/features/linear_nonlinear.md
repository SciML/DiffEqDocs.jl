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
linsolve!(x,A,b,matrix_updated=false)
```

This is an in-place function which updates `x` by solving `Ax=b`. `matrix_updated`
determines whether the matrix `A` has changed from the last call. This can be
used to smartly cache factorizations.

## Basic linsolve method: Factorization

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

## How LinSolveFactorize Was Created

In order to make your own `linsolve` functions, let's look at how the `LinSolveFactorize`
function is created. For example, for an LU-Factorization, we would like to use
`lufact!` to do our linear solving. We can directly write this as:

```julia
function linsolve!(x,A,b,update_matrix=false)
  _A = lufact!(A)
  A_ldiv_B!(x,_A,b)
end
```

This method works fine and you can pass it to the methods like

```julia
Rosenbrock23(linsolve=linsolve!)
```

and it will work, but this method does not cache `_A`, the factorization. This
means that, even if `A` has not changed, it will re-factorize the matrix.

To change this, we can instead create a call-overloaded type. The generalized form
of this is:

```julia
type LinSolveFactorize{F}
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
linsolve = LinSolveFactorize(lufact!)
```

`LinSolveFactorize` is a type which holds the factorization method and the pre-factorized
matrix. If `matrix_updated` is true, it will re-compute the factorization. Otherwise
it just solves the linear system with the cached factorization. This general
idea of using a call-overloaded type can be employed to do many other things.

## Nonlinear Solvers

Choice of nonlinear solvers is currently not available. For the current state of this work,
see [this Discourse thread](https://discourse.julialang.org/t/a-unified-interface-for-rootfinding/698/16).
