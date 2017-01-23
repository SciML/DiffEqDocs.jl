# Specifying (Non)Linear Solvers

One of the key features of DifferentialEquations.jl is its flexibility. Keeping
with this trend, many of the native Julia solvers provided by DifferentialEquations.jl
allow you to choose the method for linear and nonlinear solving. This section
details how to make that choice.

## Linear Solvers: `factorization`

For differential equation integrators which use linear solvers, an argument
to the method `factorization` determines the factorization object which is used.
For example, the `Rosenbrock23` takes in a factorization function, which we
can choose to be a QR-factorization by:

```julia
Rosenbrock23(factorization=qrfact!)
```

`factorization` is a function which returns an object that can `\`. Direct methods
like `qrfact!` will automatically cache the factorization, making it efficient
for small dense problems.

However, for large sparse problems, you can let `\` be an iterative method. For
example, using PETSc.jl, we can define our factorization function to be:

```julia
factorization = (A) -> KSP(A, ksp_type="gmres", ksp_rtol=1e-6)
```

This function creates a `KSP` type which makes `\` perform the GMRES iterative
method provided by PETSc.jl. Thus if we pass this function into the algorithm
as the factorization method, all internal linear solves will happen by PETSc.jl.

To make this compatible with more linear solvers, all that needs to happen is
such a type needs to be implemented. Feel free to contact me if you would like
help submitting a PR to help more packages use this interface.

## Nonlinear Solvers

Choice of nonlinear solvers is currently not available. For the current state of this work,
see [this Discourse thread](https://discourse.julialang.org/t/a-unified-interface-for-rootfinding/698/16).
