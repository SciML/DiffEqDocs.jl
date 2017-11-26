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
linsolve! = linsolve(Val{:init},f,x)
linsolve!(x,A,b,matrix_updated=false)
```

This is an in-place function which updates `x` by solving `Ax=b`. The user should
specify the function `linsolve(Val{:init},f,x)` which returns a `linsolve!` function.
The setting `matrix_updated` determines whether the matrix `A` has changed from the 
last call. This can be used to smartly cache factorizations.

### Basic linsolve method: Factorization

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

### How LinSolveFactorize Was Created

In order to make your own `linsolve` functions, let's look at how the `LinSolveFactorize`
function is created. For example, for an LU-Factorization, we would like to use
`lufact!` to do our linear solving. We can directly write this as:

```julia
function linsolve!(::Type{Val{:init}},f,u0)
  function _linsolve!(x,A,b,update_matrix=false)
    _A = lufact!(A)
    A_ldiv_B!(x,_A,b)
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

## Nonlinear Solvers: `nlsolve` Specification

Nonlinear solvers can be chosen via the `nlsolve` option. An `nlsolve` function
should have two dispatches:

- `nlsolve(Val{init},f,u0_prototype)` : Does an initialization phase. Returns a
  type `init_f` for later use in the solver. `u0_prototype` is the expected type
  for the initial condition `u0`.
- `nlsolve(init_f,u0)` : Solves for the root units the initialized `f` and the initial
  condition `u0`. Returns the zeros of the equation.

### Basic nlsolve mehtod: `NLSOLVEJL_SETUP`

By default, a basic nonlinear solver setup is given as `NLSOLVEJL_SETUP`. For example,
the default `nlsolve` in `Trapezoid` is

```julia
Trapezoid(nlsolve=NLSOLVEJL_SETUP())
```

This will use NLsolve.jl with autodifferentiation to solve the nonlinear systems.
`NLSOLVEJL_SETUP` has two options:

- `chunk_size` : The autodifferentiation chunk size. Integer. Defaults to ForwardDiff.jl's
  auto-detection.
- `autodiff` : Whether to use autodifferentiation. Defaults to true.

For example, to turn off autodifferentiation, use

```julia
Trapezoid(nlsolve=NLSOLVEJL_SETUP(autodiff=false))
```

### How NLSOLVEJL_SETUP Was Created

To create a nonlinear solver, you need to define the two functions. Here we use
a call-overloaded type so that way we can hold the chunk size and autodifferentiation
information.

```julia
immutable NLSOLVEJL_SETUP{CS,AD} end
Base.@pure NLSOLVEJL_SETUP(;chunk_size=0,autodiff=true) = NLSOLVEJL_SETUP{chunk_size,autodiff}()
```

The solver function just calls NLsolve and returns the zeros

```julia
(p::NLSOLVEJL_SETUP)(f,u0) = (res=NLsolve.nlsolve(f,u0); res.zero)
```

while the initialization function has a different initialization for autodifferentiation
or not:

```julia
function (p::NLSOLVEJL_SETUP{CS,AD}){CS,AD}(::Type{Val{:init}},f,u0_prototype)
  if AD
    return non_autodiff_setup(f,u0_prototype)
  else
    return autodiff_setup(f,u0_prototype,Val{determine_chunksize(initial_x,CS)})
  end
end
```

We need to declare the `get_chunksize` trait for the solver:

```julia
get_chunksize{CS,AD}(x::NLSOLVEJL_SETUP{CS,AD}) = CS
```

The initialization functions are directly for NLsolve. See the NLsolve.jl docs
for the types of inputs it expects to see. This does exactly that:

```julia
function autodiff_setup{CS}(f!, initial_x::Vector,chunk_size::Type{Val{CS}})

    permf! = (fx, x) -> f!(x, fx)

    fx2 = copy(initial_x)
    jac_cfg = ForwardDiff.JacobianConfig{CS}(initial_x, initial_x)
    g! = (x, gx) -> ForwardDiff.jacobian!(gx, permf!, fx2, x, jac_cfg)

    fg! = (x, fx, gx) -> begin
        jac_res = DiffBase.DiffResult(fx, gx)
        ForwardDiff.jacobian!(jac_res, permf!, fx2, x, jac_cfg)
        DiffBase.value(jac_res)
    end

    return DifferentiableMultivariateFunction(f!, g!, fg!)
end

function non_autodiff_setup(f!, initial_x::Vector)
  DifferentiableMultivariateFunction(f!)
end
```
