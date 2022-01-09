# [Solving Stiff Equations](@id stiff)

This tutorial is for getting into the extra features for solving stiff ordinary
differential equations in an efficient manner. Solving stiff ordinary
differential equations requires specializing the linear solver on properties of
the Jacobian in order to cut down on the ``\mathcal{O}(n^3)`` linear solve and
the ``\mathcal{O}(n^2)`` back-solves. Note that these same functions and
controls also extend to stiff SDEs, DDEs, DAEs, etc.

!!! note

    This tutorial is for advanced users to dive into advanced features!
    DifferentialEquations.jl automates most of this usage, so we recommend
    users try `solve(prob)` with the automatic algorithm first!

## Code Optimization for Differential Equations

Solving stiff differential equations requires speed. Here's a few things to
keep in mind.

### Writing Efficient Code

For a detailed tutorial on how to optimize one's DifferentialEquations.jl code,
please see the
[Optimizing DiffEq Code tutorial](http://tutorials.juliadiffeq.org/html/introduction/03-optimizing_diffeq_code.html).

### Choosing a Good Solver

Choosing a good solver is required for getting top notch speed. General
recommendations can be found on the solver page (for example, the
[ODE Solver Recommendations](@ref ode_solve)).
The current recommendations can be simplified to a Rosenbrock method
(`Rosenbrock23` or `Rodas5`) for smaller (<50 ODEs) problems, ESDIRK methods
for slightly larger (`TRBDF2` or `KenCarp4` for <2000 ODEs), and `QNDF` for even
larger problems. `lsoda` from [LSODA.jl](https://github.com/rveltz/LSODA.jl) is
sometimes worth a try for the medium sized category.

More details on the solver to choose can be found by benchmarking. See the
[SciMLBenchmarks](https://github.com/JuliaDiffEq/SciMLBenchmarks.jl) to
compare many solvers on many problems.

### Check Out the Speed FAQ

See [this FAQ](@ref faq_performance)
for information on common pitfalls and how to improve performance.

### Setting Up Your Julia Installation for Speed

Julia uses an underlying BLAS implementation for its matrix multiplications
and factorizations. This library is automatically multithreaded and accelerates
the internal linear algebra of DifferentialEquations.jl. However, for optimality,
you should make sure that the number of BLAS threads that you are using matches
the number of physical cores and not the number of logical cores. See
[this issue for more details](https://github.com/JuliaLang/julia/issues/33409).

To check the number of BLAS threads, use:

```julia
ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
```

If I want to set this directly to 4 threads, I would use:

```julia
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(4)
```

Additionally, in some cases Intel's MKL might be a faster BLAS than the standard
BLAS that ships with Julia (OpenBLAS). To switch your BLAS implementation, you
can use [MKL.jl](https://github.com/JuliaComputing/MKL.jl) which will accelerate
the linear algebra routines. This is done via:

```julia
using MKL
```

### Use Accelerator Hardware

When possible, use GPUs. If your ODE system is small and you need to solve it
with very many different parameters, see the
[ensembles interface](@ref ensemble)
and [DiffEqGPU.jl](https://github.com/JuliaDiffEq/DiffEqGPU.jl). If your problem
is large, consider using a [CuArray](https://github.com/JuliaGPU/CuArrays.jl)
for the state to allow for GPU-parallelism of the internal linear algebra.

## Speeding Up Jacobian Calculations

When one is using an implicit or semi-implicit differential equation solver,
the Jacobian must be built at many iterations and this can be one of the most
expensive steps. There are two pieces that must be optimized in order to reach
maximal efficiency when solving stiff equations: the sparsity pattern and the
construction of the Jacobian. The construction is filling the matrix
`J` with values, while the sparsity pattern is what `J` to use.

The sparsity pattern is given by a prototype matrix, the `jac_prototype`, which
will be copied to be used as `J`. The default is for `J` to be a `Matrix`,
i.e. a dense matrix. However, if you know the sparsity of your problem, then
you can pass a different matrix type. For example, a `SparseMatrixCSC` will
give a sparse matrix. Additionally, structured matrix types like `Tridiagonal`,
`BandedMatrix` (from
[BandedMatrices.jl](https://github.com/JuliaMatrices/BandedMatrices.jl)),
`BlockBandedMatrix` (from
[BlockBandedMatrices.jl](https://github.com/JuliaMatrices/BlockBandedMatrices.jl)),
and more can be given. DifferentialEquations.jl will internally use this matrix
type, making the factorizations faster by utilizing the specialized forms.

For the construction, there are 3 ways to fill `J`:

- The default, which uses normal finite/automatic differentiation
- A function `jac(J,u,p,t)` which directly computes the values of `J`
- A `colorvec` which defines a sparse differentiation scheme.

We will now showcase how to make use of this functionality with growing complexity.

### Declaring Jacobian Functions

Let's solve the Rober equations:

```math
\begin{aligned}
\frac{dy_1}{dt} &= -0.04y₁ + 10^4 y_2 y_3 \\
\frac{dy_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
\frac{dy_3}{dt} &= 3*10^7 y_{2}^2 \\
\end{aligned}
```

In order to reduce the Jacobian construction cost, one can describe a Jacobian
function by using the `jac` argument for the `ODEFunction`. First, let's do
a standard `ODEProblem`:

```julia
using DifferentialEquations
function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁+k₃*y₂*y₃
  du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du[3] =  k₂*y₂^2
  nothing
end
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),[0.04,3e7,1e4])
sol = solve(prob)
plot(sol,tspan=(1e-2,1e5),xscale=:log10)
```

![IntroDAEPlot](../assets/intro_dae_plot.png)

```julia-repl
julia> using BenchmarkTools
julia> @btime solve(prob)
97.000 μs (1832 allocations: 132.30 KiB)
```

Now we want to add the Jacobian. First we have to derive the Jacobian
``\frac{df_i}{du_j}`` which is `J[i,j]`. From this we get:

```julia
function rober_jac(J,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  J[1,1] = k₁ * -1
  J[2,1] = k₁
  J[3,1] = 0
  J[1,2] = y₃ * k₃
  J[2,2] = y₂ * k₂ * -2 + y₃ * k₃ * -1
  J[3,2] = y₂ * 2 * k₂
  J[1,3] = k₃ * y₂
  J[2,3] = k₃ * y₂ * -1
  J[3,3] = 0
  nothing
end
f = ODEFunction(rober, jac=rober_jac)
prob_jac = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
```
```julia-repl
julia> @btime solve(prob_jac)
74.500 μs (1475 allocations: 100.50 KiB)
```

### Automatic Derivation of Jacobian Functions

But that was hard! If you want to take the symbolic Jacobian of numerical
code, we can make use of [ModelingToolkit.jl](https://github.com/JuliaDiffEq/ModelingToolkit.jl)
to symbolicify the numerical code and do the symbolic calculation and return
the Julia code for this.

```julia
using ModelingToolkit
de = modelingtoolkitize(prob)
ModelingToolkit.generate_jacobian(de)[2] # Second is in-place
```

which outputs:

```julia
julia> ModelingToolkit.generate_jacobian(de)[2] # Second is in-place
:(function (ˍ₋out, ˍ₋arg1, ˍ₋arg2, t)
      #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:303 =#
      #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:304 =#
      let var"x₁(t)" = #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:190 =# @inbounds(ˍ₋arg1[1]), var"x₂(t)" = #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:190 =# @inbounds(ˍ₋arg1[2]), var"x₃(t)" = #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:190 =# @inbounds(ˍ₋arg1[3]), α₁ = #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:190 =# @inbounds(ˍ₋arg2[1]), α₂ = #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:190 =# @inbounds(ˍ₋arg2[2]), α₃ = #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:190 =# @inbounds(ˍ₋arg2[3])
          #= C:\Users\accou\.julia\dev\Symbolics\src\build_function.jl:378 =#
          #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:350 =# @inbounds begin
                  #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:346 =#
                  ˍ₋out[1] = (*)(-1, α₁)
                  ˍ₋out[2] = α₁
                  ˍ₋out[3] = 0
                  ˍ₋out[4] = (*)(α₃, var"x₃(t)")
                  ˍ₋out[5] = (+)((*)((*)(-2, α₂), var"x₂(t)"), (*)((*)(-1, α₃), var"x₃(t)"))
                  ˍ₋out[6] = (*)((*)(2, α₂), var"x₂(t)")
                  ˍ₋out[7] = (*)(α₃, var"x₂(t)")
                  ˍ₋out[8] = (*)((*)(-1, α₃), var"x₂(t)")
                  ˍ₋out[9] = 0
                  #= C:\Users\accou\.julia\packages\SymbolicUtils\0KTj4\src\code.jl:348 =#
                  nothing
              end
      end
  end)
```

Now let's use that to give the analytical solution Jacobian:

```julia
jac = eval(ModelingToolkit.generate_jacobian(de)[2])
f = ODEFunction(rober, jac=jac)
prob_jac = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
```

### Accelerating Small ODE Solves with Static Arrays

If the ODE is sufficiently small (<20 ODEs or so), using [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)
for the state variables can greatly enhance the performance. This is done by
making `u0` a `StaticArray` and writing an out-of-place non-mutating dispatch
for static arrays, for the ROBER problem, this looks like:

```julia
using DifferentialEquations, StaticArrays
function rober_static(u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du1 = -k₁*y₁+k₃*y₂*y₃
  du2 =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du3 =  k₂*y₂^2
  SA[du1,du2,du3]
end
prob = ODEProblem(rober_static,SA[1.0,0.0,0.0],(0.0,1e5),SA[0.04,3e7,1e4])
sol = solve(prob,Rosenbrock23())
```

If we benchmark this we see a really fast solution with really low allocation
counts:

```julia
@btime sol = solve(prob,Rosenbrock23())
# 12.100 μs (87 allocations: 18.53 KiB)
```

This version is thus very amenable to multithreading and other forms of parallelism.

## Accelerating Large Stiff ODE / PDE Solves

This section is focused on the tools and tricks for getting the most performance
out of very large ODE solves, such as those resulting from PDE discretizations.

### Declaring a Sparse Jacobian

Jacobian sparsity is declared by the `jac_prototype` argument in the `ODEFunction`.
Note that you should only do this if the sparsity is high, for example, 0.1%
of the matrix is non-zeros, otherwise the overhead of sparse matrices can be higher
than the gains from sparse differentiation!

But as a demonstration, let's build a sparse matrix for the Rober problem. We
can do this by gathering the `I` and `J` pairs for the non-zero components, like:

```julia
I = [1,2,1,2,3,1,2]
J = [1,1,2,2,2,3,3]

using SparseArrays
jac_prototype = sparse(I,J,1.0)
```

Now this is the sparse matrix prototype that we want to use in our solver, which
we then pass like:

```julia
f = ODEFunction(rober, jac=jac, jac_prototype=jac_prototype)
prob_jac = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
```

### Automatic Sparsity Detection

One of the useful companion tools for DifferentialEquations.jl is
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).
This allows for automatic declaration of Jacobian sparsity types. To see this
in action, let's look at the 2-dimensional Brusselator equation:

```julia
const N = 32
const xyd_brusselator = range(0,stop=1,length=N)
brusselator_f(x, y, t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.
limit(a, N) = a == N+1 ? 1 : a == 0 ? N : a
function brusselator_2d_loop(du, u, p, t)
  A, B, alpha, dx = p
  alpha = alpha/dx^2
  @inbounds for I in CartesianIndices((N, N))
    i, j = Tuple(I)
    x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
    ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
    du[i,j,1] = alpha*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) +
                B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1] + brusselator_f(x, y, t)
    du[i,j,2] = alpha*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
    end
end
p = (3.4, 1., 10., step(xyd_brusselator))
```

Given this setup, we can give and example `du` and `u` and call `jacobian_sparsity`
on our function with the example arguments and it will kick out a sparse matrix
with our pattern, that we can turn into our `jac_prototype`.

```julia
using Symbolics
du0 = copy(u0)
jac_sparsity = Symbolics.jacobian_sparsity((du,u)->brusselator_2d_loop(du,u,p,0.0),du0,u0)

2048×2048 SparseArrays.SparseMatrixCSC{Bool, Int64} with 12288 stored entries:
⠻⣦⡀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀
⠀⠈⠻⣦⡀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀
⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠈⠳⣄⠀⠀
⡀⠀⠀⠀⠀⠈⠻⣦⠀⠀⠀⠀⠀⠈⠳⣄
⠙⢦⡀⠀⠀⠀⠀⠀⠻⣦⡀⠀⠀⠀⠀⠈
⠀⠀⠙⢦⡀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀
⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠈⠻⣦⡀⠀
⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠈⠻⣦
```

Notice Julia gives a nice print out of the sparsity pattern. That's neat, and
would be tedious to build by hand! Now we just pass it to the `ODEFunction`
like as before:

```julia
f = ODEFunction(brusselator_2d_loop;jac_prototype=float.(jac_sparsity))
```

Build the `ODEProblem`:

```julia
function init_brusselator_2d(xyd)
  N = length(xyd)
  u = zeros(N, N, 2)
  for I in CartesianIndices((N, N))
    x = xyd[I[1]]
    y = xyd[I[2]]
    u[I,1] = 22*(y*(1-y))^(3/2)
    u[I,2] = 27*(x*(1-x))^(3/2)
  end
  u
end
u0 = init_brusselator_2d(xyd_brusselator)
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop,
                                     u0,(0.,11.5),p)

prob_ode_brusselator_2d_sparse = ODEProblem(f,
                                     u0,(0.,11.5),p)
```

Now let's see how the version with sparsity compares to the version without:

```julia
@btime solve(prob_ode_brusselator_2d,TRBDF2(),save_everystep=false) # 2.771 s (5452 allocations: 65.73 MiB)
@btime solve(prob_ode_brusselator_2d_sparse,TRBDF2(),save_everystep=false) # 680.612 ms (37905 allocations: 359.34 MiB)
```

## Defining Linear Solver Routines and Jacobian-Free Newton-Krylov

A completely different way to optimize the linear solvers for large sparse
matrices is to use a Krylov subpsace method. This requires choosing a linear
solver for changing to a Krylov method. Optionally, one can use a Jacobian-free
operator to reduce the memory requirements.

### Declaring a Jacobian-Free Newton-Krylov Implementation

To swap the linear solver out, we use the `linsolve` command and choose the
GMRES linear solver.

```julia
@btime solve(prob_ode_brusselator_2d,KenCarp47(linsolve=KrylovJL_GMRES()),save_everystep=false)
# 707.439 ms (173868 allocations: 31.07 MiB)
```

Notice that this acceleration does not require the definition of a sparsity
pattern and can thus be an easier way to scale for large problems. For more
information on linear solver choices, see the
[linear solver documentation](@ref linear_nonlinear). `linsolve` choices are any
valid [LinearSolve.jl](http://linearsolve.sciml.ai/dev/) solver.

### Adding a Preconditioner

Any [LinearSolve.jl-compatible preconditioner](http://linearsolve.sciml.ai/dev/basics/Preconditioners/)
can be used as a preconditioner in the linear solver interface. To define
preconditioners, one must define a `precs` function in compatible stiff ODE
solvers which returns the left and right preconditioners, matrices which
approximate the inverse of `W = I - gamma*J` used in the solution of the ODE.
An example of this with using [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl)
is as follows:

```julia
using IncompleteLU
function incompletelu(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = ilu(convert(AbstractMatrix,W), τ = 50.0)
  else
    Pl = Plprev
  end
  Pl,nothing
end

# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv,Ti}) where {Tv,Ti} = Tv

@time solve(prob_ode_brusselator_2d_sparse,KenCarp47(linsolve=KrylovJL_GMRES(),precs=incompletelu,concrete_jac=true),save_everystep=false);
# 174.386 ms (61756 allocations: 61.38 MiB)
```

Notice a few things about this preconditioner. This preconditioner uses the
sparse Jacobian, and thus we set `concrete_jac=true` to tell the algorithm to
generate the Jacobian (otherwise, a Jacobian-free algorithm is used with GMRES
by default). Then `newW = true` whenever a new `W` matrix is computed, and
`newW=nothing` during the startup phase of the solver. Thus we do a check
`newW === nothing || newW` and when true, it's only at these points when we
we update the preconditioner, otherwise we just pass on the previous version.
We use `convert(AbstractMatrix,W)` to get the concrete `W` matrix (matching
`jac_prototype`, thus `SpraseMatrixCSC`) which we can use in the preconditioner's
definition. Then we use `IncompleteLU.ilu` on that sparse matrix to generate
the preconditioner. We return `Pl,nothing` to say that our preconditioner is a
left preconditioner, and that there is no right preconditioning.

This method thus uses both the Krylov solver and the sparse Jacobian. Not only
that, it is faster than both implementations! IncompleteLU is fussy in that it
requires a well-tuned `τ` parameter. Another option is to use
[AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
which is more automatic. The setup is very similar to before:

```julia
using AlgebraicMultigrid
function algebraicmultigrid(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = aspreconditioner(ruge_stuben(convert(AbstractMatrix,W)))
  else
    Pl = Plprev
  end
  Pl,nothing
end

# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::AlgebraicMultigrid.Preconditioner) = Float64

@btime solve(prob_ode_brusselator_2d_sparse,KenCarp47(linsolve=KrylovJL_GMRES(),precs=algebraicmultigrid,concrete_jac=true),save_everystep=false);
# 372.528 ms (61179 allocations: 160.82 MiB)
```

or with a Jacobi smoother:

```julia
function algebraicmultigrid2(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    A = convert(AbstractMatrix,W)
    Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(A, presmoother = AlgebraicMultigrid.Jacobi(rand(size(A,1))), postsmoother = AlgebraicMultigrid.Jacobi(rand(size(A,1)))))
  else
    Pl = Plprev
  end
  Pl,nothing
end

@btime solve(prob_ode_brusselator_2d_sparse,KenCarp47(linsolve=KrylovJL_GMRES(),precs=algebraicmultigrid2,concrete_jac=true),save_everystep=false);
# 293.476 ms (65714 allocations: 170.23 MiB)
```

For more information on the preconditioner interface, see the
[linear solver documentation](@ref linear_nonlinear).

## Using Structured Matrix Types

If your sparsity pattern follows a specific structure, for example a banded
matrix, then you can declare `jac_prototype` to be of that structure and then
additional optimizations will come for free. Note that in this case, it is
not necessary to provide a `colorvec` since the color vector will be analytically
derived from the structure of the matrix.

The matrices which are allowed are those which satisfy the
[ArrayInterface.jl](https://github.com/JuliaDiffEq/ArrayInterface.jl) interface
for automatically-colorable matrices. These include:

- Bidiagonal
- Tridiagonal
- SymTridiagonal
- BandedMatrix ([BandedMatrices.jl](https://github.com/JuliaMatrices/BandedMatrices.jl))
- BlockBandedMatrix ([BlockBandedMatrices.jl](https://github.com/JuliaMatrices/BlockBandedMatrices.jl))

Matrices which do not satisfy this interface can still be used, but the matrix
coloring will not be automatic, and an appropriate linear solver may need to
be given (otherwise it will default to attempting an LU-decomposition).

## Sundials-Specific Handling

While much of the setup makes the transition to using Sundials automatic, there
are some differences between the pure Julia implementations and the Sundials
implementations which must be taken note of. These are all detailed in the
[Sundials solver documentation](http://docs.juliadiffeq.org/dev/solvers/ode_solve#Sundials.jl-1),
but here we will highlight the main details which one should make note of.

Defining a sparse matrix and a Jacobian for Sundials works just like any other
package. The core difference is in the choice of the linear solver. With Sundials,
the linear solver choice is done with a Symbol in the `linear_solver` from a
preset list. Particular choices of note are `:Band` for a banded matrix and
`:GMRES` for using GMRES. If you are using Sundials, `:GMRES` will not require
defining the JacVecOperator, and instead will always make use of a Jacobian-Free
Newton Krylov (with numerical differentiation). Thus on this problem we could do:

```julia
using Sundials
@btime solve(prob_ode_brusselator_2d,CVODE_BDF(),save_everystep=false) # 13.280 s (51457 allocations: 2.43 MiB)
# Simplest speedup: use :LapackDense
@btime solve(prob_ode_brusselator_2d,CVODE_BDF(linear_solver=:LapackDense),save_everystep=false) # 2.024 s (51457 allocations: 2.43 MiB)
# GMRES Version: Doesn't require any extra stuff!
@btime solve(prob_ode_brusselator_2d,CVODE_BDF(linear_solver=:GMRES),save_everystep=false) # 213.800 ms (58353 allocations: 2.64 MiB)
```

Notice that using sparse matrices with Sundials requires an analytical Jacobian
function. We will use [ModelingToolkit.jl](https://mtk.sciml.ai/dev/)'s
`modelingtoolkitize` to automatically generate this:

```julia
using ModelingToolkit
prob_ode_brusselator_2d_mtk = ODEProblem(modelingtoolkitize(prob_ode_brusselator_2d_sparse),[],(0.0,11.5),jac=true,sparse=true);
@btime solve(prob_ode_brusselator_2d_mtk,CVODE_BDF(linear_solver=:KLU),save_everystep=false) # 493.908 ms (1358 allocations: 5.79 MiB)
```

### Using Preconditioners with Sundials

Details for setting up a preconditioner with Sundials can be found at the
[Sundials solver page](@ref ode_solve_sundials). Sundials algorithms are very
different from the standard Julia-based algorithms in that they require the
user does all handling of the Jacobian matrix. To do this, you must define a
`psetup` function that sets up the preconditioner and then a `prec` function
that is the action of the preconditioner on a vector. For the `psetup` function,
we need to first compute the `W = I - gamma*J` matrix before computing the
preconditioner on it. For the ILU example above, this is done for Sundials like:

```julia
using LinearAlgebra
u0 = prob_ode_brusselator_2d_mtk.u0
p  = prob_ode_brusselator_2d_mtk.p
const jaccache = prob_ode_brusselator_2d_mtk.f.jac(u0,p,0.0)
const W = I - 1.0*jaccache

prectmp = ilu(W, τ = 50.0)
const preccache = Ref(prectmp)

function psetupilu(p, t, u, du, jok, jcurPtr, gamma)
  if jok
    prob_ode_brusselator_2d_mtk.f.jac(jaccache,u,p,t)
    jcurPtr[] = true

    # W = I - gamma*J
    @. W = -gamma*jaccache
    idxs = diagind(W)
    @. @view(W[idxs]) = @view(W[idxs]) + 1

    # Build preconditioner on W
    preccache[] = ilu(W, τ = 5.0)
  end
end
```

Then the preconditioner action is to simply use the `ldiv!` of the generated
preconditioner:

```julia
function precilu(z,r,p,t,y,fy,gamma,delta,lr)
  ldiv!(z,preccache[],r)
end
```

We then simply pass these functions to the Sundials solver with a choice of
`prec_side=1` to indicate that it is a left-preconditioner:

```julia
@btime solve(prob_ode_brusselator_2d_sparse,CVODE_BDF(linear_solver=:GMRES,prec=precilu,psetup=psetupilu,prec_side=1),save_everystep=false);
# 87.176 ms (17717 allocations: 77.08 MiB)
```

And similarly for algebraic multigrid:

```julia
prectmp2 = aspreconditioner(ruge_stuben(W, presmoother = AlgebraicMultigrid.Jacobi(rand(size(W,1))), postsmoother = AlgebraicMultigrid.Jacobi(rand(size(W,1)))))
const preccache2 = Ref(prectmp2)
function psetupamg(p, t, u, du, jok, jcurPtr, gamma)
  if jok
    prob_ode_brusselator_2d_mtk.f.jac(jaccache,u,p,t)
    jcurPtr[] = true

    # W = I - gamma*J
    @. W = -gamma*jaccache
    idxs = diagind(W)
    @. @view(W[idxs]) = @view(W[idxs]) + 1

    # Build preconditioner on W
    preccache2[] = aspreconditioner(ruge_stuben(W, presmoother = AlgebraicMultigrid.Jacobi(rand(size(W,1))), postsmoother = AlgebraicMultigrid.Jacobi(rand(size(W,1)))))
  end
end

function precamg(z,r,p,t,y,fy,gamma,delta,lr)
  ldiv!(z,preccache2[],r)
end

@btime solve(prob_ode_brusselator_2d_sparse,CVODE_BDF(linear_solver=:GMRES,prec=precamg,psetup=psetupamg,prec_side=1),save_everystep=false);
# 136.431 ms (30682 allocations: 275.68 MiB)
```

## Handling Mass Matrices

Instead of just defining an ODE as ``u' = f(u,p,t)``, it can be common to express
the differential equation in the form with a mass matrix:

```math
Mu' = f(u,p,t)
```

where ``M`` is known as the mass matrix. Let's solve the Robertson equation.
At the top we wrote this equation as:

```math
\begin{aligned}
dy_1 &= -0.04 y_1 + 10^4 y_2 y_3 \\
dy_2 &=  0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
dy_3 &= 3*10^7 y_{2}^2 \\
\end{aligned}
```

But we can instead write this with a conservation relation:

```math
\begin{aligned}
\frac{dy_1}{dt} &= -0.04 y_1 + 10^4 y_2 y_3 \\
\frac{dy_2}{dt} &=  0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
1 &=  y_{1} + y_{2} + y_{3} \\
\end{aligned}
```

In this form, we can write this as a mass matrix ODE where ``M`` is singular
(this is another form of a differential-algebraic equation (DAE)). Here, the
last row of `M` is just zero. We can implement this form as:

```julia
using DifferentialEquations
function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁ + k₃*y₂*y₃
  du[2] =  k₁*y₁ - k₃*y₂*y₃ - k₂*y₂^2
  du[3] =  y₁ + y₂ + y₃ - 1
  nothing
end
M = [1. 0  0
     0  1. 0
     0  0  0]
f = ODEFunction(rober,mass_matrix=M)
prob_mm = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol = solve(prob_mm,Rodas5(),reltol=1e-8,abstol=1e-8)

plot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))
```

![IntroDAEPlot](../assets/intro_dae_plot.png)

!!! note
    If your mass matrix is singular, i.e. your system is a DAE, then you
    need to make sure you choose
    [a solver that is compatible with DAEs](@ref dae_solve_full)
