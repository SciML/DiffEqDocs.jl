# [Solving Large Stiff Equations](@id stiff)

This tutorial is for getting into the extra features for solving large stiff ordinary
differential equations efficiently. Solving stiff ordinary
differential equations requires specializing the linear solver on properties of
the Jacobian in order to cut down on the ``\mathcal{O}(n^3)`` linear solve and
the ``\mathcal{O}(n^2)`` back-solves. Note that these same functions and
controls also extend to stiff SDEs, DDEs, DAEs, etc. This tutorial is for large-scale
models, such as those derived for semi-discretizations of partial differential
equations (PDEs). For example, we will use the stiff Brusselator partial
differential equation (BRUSS).

!!! note
    
    This tutorial is for advanced users to dive into advanced features!
    DifferentialEquations.jl automates most of this usage, so we recommend
    users try `solve(prob)` with the automatic algorithm first!

## Definition of the Brusselator Equation

!!! note
    
    Feel free to skip this section: it simply defines the example problem.

The Brusselator PDE is defined on a unit square periodic domain as follows:

```math
\begin{align}
\frac{\partial U}{\partial t} &= 1 + U^2V - 4.4U + \alpha \nabla^2 U + f(x, y, t),\\
\frac{\partial V}{\partial t} &= 3.4U - U^2V + \alpha \nabla^2 V,
\end{align}
```

where

```math
f(x, y, t) = \begin{cases}
5 & \quad \text{if } (x-0.3)^2+(y-0.6)^2 ≤ 0.1^2 \text{ and } t ≥ 1.1\\
0 & \quad \text{else}
\end{cases}, \mathrm{and}
```

$\nabla^2 = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$ is the two dimensional Laplacian operator. The above equations are to be solved for a time interval $t \in [0, 11.5]$ subject to the initial conditions

```math
\begin{align}
U(x, y, 0) &= 22\cdot (y(1-y))^{3/2} \\
V(x, y, 0) &= 27\cdot (x(1-x))^{3/2}
\end{align},
```

and the periodic boundary conditions

```math
\begin{align}
U(x+1,y,t) &= U(x,y,t) \\
V(x,y+1,t) &= V(x,y,t).
\end{align}
```

To solve this PDE, we will discretize it into a system of ODEs with the finite
difference method. We discretize the unit square domain with `N` grid points in each direction.
`U[i,j]` and `V[i,j]` then represent the value of the discretized field at a given point in time, i.e.

```
U[i,j] = U(i*dx,j*dy)
V[i,j] = V(i*dx,j*dy)
```

where `dx = dy = 1/N`. To implement our ODE system, we collect both `U` and `V` in a single array `u` of size `(N,N,2)` with `u[i,j,1] = U[i,j]` and `u[i,j,2] = V[i,j]`. This approach can be easily generalized to PDEs with larger number of field variables.

Using a three-point stencil, the Laplacian operator discretizes into a tridiagonal matrix with elements `[1 -2 1]` and a `1` in the top, bottom, left, and right corners coming from the periodic boundary conditions. The nonlinear terms are implemented pointwise in a straightforward manner.

The resulting `ODEProblem` definition is:

```@example stiff1
import DifferentialEquations as DE, LinearAlgebra, SparseArrays

const N = 32
const xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a
function brusselator_2d_loop(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
        limit(j - 1, N)
        du[i, j, 1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)
        du[i, j, 2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                       4u[i, j, 2]) +
                      A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
end
p = (3.4, 1.0, 10.0, step(xyd_brusselator))

function init_brusselator_2d(xyd)
    N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    u
end
u0 = init_brusselator_2d(xyd_brusselator)
prob_ode_brusselator_2d = DE.ODEProblem(brusselator_2d_loop, u0, (0.0, 11.5), p)
```

## Choosing Jacobian Types

When one is using an implicit or semi-implicit differential equation solver,
the Jacobian must be built at many iterations, and this can be one of the most
expensive steps. There are two pieces that must be optimized in order to reach
maximal efficiency when solving stiff equations: the sparsity pattern and the
construction of the Jacobian. The construction is filling the matrix
`J` with values, while the sparsity pattern is what `J` to use.

The sparsity pattern is given by a prototype matrix, the `jac_prototype`, which
will be copied to be used as `J`. The default is for `J` to be a `Matrix`,
i.e. a dense matrix. However, if you know the sparsity of your problem, then
you can pass a different matrix type. For example, a `SparseMatrixCSC` will
give a sparse matrix. Other sparse matrix types include:

  - Bidiagonal
  - Tridiagonal
  - SymTridiagonal
  - BandedMatrix ([BandedMatrices.jl](https://github.com/JuliaLinearAlgebra/BandedMatrices.jl))
  - BlockBandedMatrix ([BlockBandedMatrices.jl](https://github.com/JuliaLinearAlgebra/BlockBandedMatrices.jl))

DifferentialEquations.jl will internally use this matrix
type, making the factorizations faster by using the specialized forms.

## Declaring a Sparse Jacobian with Automatic Sparsity Detection

Jacobian sparsity is declared by the `jac_prototype` argument in the `ODEFunction`.
Note that you should only do this if the sparsity is high, for example, 0.1%
of the matrix is non-zeros, otherwise the overhead of sparse matrices can be higher
than the gains from sparse differentiation!

[ADTypes.jl](https://github.com/SciML/ADTypes.jl) provides a [common interface for automatic sparsity detection](https://sciml.github.io/ADTypes.jl/stable/#Sparsity-detector)
via its function `jacobian_sparsity`.
This function can be called using sparsity detectors from [SparseConnectivityTracer.jl](https://github.com/adrhill/SparseConnectivityTracer.jl)
or [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

We can give an example `du` and `u` and call `jacobian_sparsity` on our function with the example arguments,
and it will kick out a sparse matrix with our pattern, that we can turn into our `jac_prototype`.

Let's try SparseConnectivityTracer's [`TracerSparsityDetector`](https://adrianhill.de/SparseConnectivityTracer.jl/stable/user/api/#SparseConnectivityTracer.TracerSparsityDetector):

```@example stiff1
import SparseConnectivityTracer, ADTypes
detector = SparseConnectivityTracer.TracerSparsityDetector()
du0 = copy(u0)
jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> brusselator_2d_loop(du, u, p, 0.0), du0, u0, detector)
```

Using a different backend for sparsity detection just requires swapping out the detector,
e.g. for Symbolics' [`SymbolicsSparsityDetector`](https://docs.sciml.ai/Symbolics/stable/manual/sparsity_detection/#Symbolics.SymbolicsSparsityDetector).

Notice that Julia gives a nice print out of the sparsity pattern.
That's neat, and would be tedious to build by hand!
Now we just pass it to the `ODEFunction` like as before:

```@example stiff1
f = DE.ODEFunction(brusselator_2d_loop; jac_prototype = float.(jac_sparsity))
```

Build the `ODEProblem`:

```@example stiff1
prob_ode_brusselator_2d_sparse = DE.ODEProblem(f, u0, (0.0, 11.5), p)
```

Now let's see how the version with sparsity compares to the version without:

```@example stiff1
import BenchmarkTools as BT # for @btime
BT.@btime DE.solve(prob_ode_brusselator_2d, DE.TRBDF2(); save_everystep = false);
BT.@btime DE.solve(prob_ode_brusselator_2d_sparse, DE.TRBDF2(); save_everystep = false);
BT.@btime DE.solve(prob_ode_brusselator_2d_sparse, DE.KenCarp47(; linsolve = DE.KLUFactorization());
    save_everystep = false);
nothing # hide
```

Note that depending on the properties of the sparsity pattern, one may want
to try alternative linear solvers such as `DE.TRBDF2(linsolve = DE.KLUFactorization())`
or `DE.TRBDF2(linsolve = DE.UMFPACKFactorization())`.

## Using Jacobian-Free Newton-Krylov

A completely different way to optimize the linear solvers for large sparse
matrices is to use a Krylov subspace method. This requires choosing a linear
solver for changing to a Krylov method. To swap the linear solver out, we use
the `linsolve` command and choose the GMRES linear solver.

```@example stiff1
BT.@btime DE.solve(prob_ode_brusselator_2d, DE.KenCarp47(; linsolve = DE.KrylovJL_GMRES());
    save_everystep = false);
nothing # hide
```

Notice that this acceleration does not require the definition of a sparsity
pattern, and can thus be an easier way to scale for large problems. For more
information on linear solver choices, see the
[linear solver documentation](@ref linear_nonlinear). `linsolve` choices are any
valid [LinearSolve.jl](https://linearsolve.sciml.ai/dev/) solver.

!!! note
    
    Switching to a Krylov linear solver will automatically change the ODE solver
    into Jacobian-free mode, dramatically reducing the memory required. This can
    be overridden by adding `concrete_jac=true` to the algorithm.

## Adding a Preconditioner

Any [LinearSolve.jl-compatible preconditioner](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/)
can be used as a preconditioner in the linear solver interface. To define
preconditioners, one must define a `precs` function in compatible stiff ODE
solvers which returns the left and right preconditioners, matrices which
approximate the inverse of `W = I - gamma*J` used in the solution of the ODE.
An example of this with using [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl)
is as follows:

```@example stiff1
import IncompleteLU
function incompletelu(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = IncompleteLU.ilu(convert(AbstractMatrix, W), τ = 50.0)
    else
        Pl = Plprev
    end
    Pl, nothing
end

# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv, Ti}) where {Tv, Ti} = Tv

BT.@btime DE.solve(prob_ode_brusselator_2d_sparse,
    DE.KenCarp47(; linsolve = DE.KrylovJL_GMRES(), precs = incompletelu,
        concrete_jac = true); save_everystep = false);
nothing # hide
```

Notice a few things about this preconditioner. This preconditioner uses the
sparse Jacobian, and thus we set `concrete_jac=true` to tell the algorithm to
generate the Jacobian (otherwise, a Jacobian-free algorithm is used with GMRES
by default). Then `newW = true` whenever a new `W` matrix is computed, and
`newW=nothing` during the startup phase of the solver. Thus, we do a check
`newW === nothing || newW` and when true, it's only at these points when
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

```@example stiff1
import AlgebraicMultigrid
function algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(convert(AbstractMatrix, W)))
    else
        Pl = Plprev
    end
    Pl, nothing
end

BT.@btime DE.solve(prob_ode_brusselator_2d_sparse,
    DE.KenCarp47(; linsolve = DE.KrylovJL_GMRES(), precs = algebraicmultigrid,
        concrete_jac = true); save_everystep = false);
nothing # hide
```

or with a Jacobi smoother:

```@example stiff1
function algebraicmultigrid2(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        A = convert(AbstractMatrix, W)
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(A,
            presmoother = AlgebraicMultigrid.Jacobi(rand(size(A,
                1))),
            postsmoother = AlgebraicMultigrid.Jacobi(rand(size(A,
                1)))))
    else
        Pl = Plprev
    end
    Pl, nothing
end

BT.@btime DE.solve(prob_ode_brusselator_2d_sparse,
    DE.KenCarp47(; linsolve = DE.KrylovJL_GMRES(), precs = algebraicmultigrid2,
        concrete_jac = true); save_everystep = false);
nothing # hide
```

For more information on the preconditioner interface, see the
[linear solver documentation](@ref linear_nonlinear).

## Sundials-Specific Handling

While much of the setup makes the transition to using Sundials automatic, there
are some differences between the pure Julia implementations and the Sundials
implementations which must be taken note of. These are all detailed in the
[Sundials solver documentation](https://docs.sciml.ai/DiffEqDocs/stable/api/sundials/),
but here we will highlight the main details which one should make note of.

Defining a sparse matrix and a Jacobian for Sundials works just like any other
package. The core difference is in the choice of the linear solver. With Sundials,
the linear solver choice is done with a Symbol in the `linear_solver` from a
preset list. Particular choices of note are `:Band` for a banded matrix and
`:GMRES` for using GMRES. If you are using Sundials, `:GMRES` will not require
defining the JacVecOperator, and instead will always make use of a Jacobian-Free
Newton Krylov (with numerical differentiation). Thus, on this problem we could do:

```@example stiff1
import Sundials
BT.@btime DE.solve(prob_ode_brusselator_2d, Sundials.CVODE_BDF(); save_everystep = false);
# Simplest speedup: use :LapackDense
BT.@btime DE.solve(prob_ode_brusselator_2d, Sundials.CVODE_BDF(; linear_solver = :LapackDense);
    save_everystep = false);
# GMRES Version: Doesn't require any extra stuff!
BT.@btime DE.solve(prob_ode_brusselator_2d, Sundials.CVODE_BDF(; linear_solver = :GMRES);
    save_everystep = false);
nothing # hide
```

Notice that using sparse matrices with Sundials requires an analytical Jacobian
function. We will use [ModelingToolkit.jl](https://mtk.sciml.ai/dev/)'s
`modelingtoolkitize` to automatically generate this:

```@example stiff1
import ModelingToolkit as MTK
prob_ode_brusselator_2d_mtk = DE.ODEProblem(
    MTK.complete(MTK.modelingtoolkitize(prob_ode_brusselator_2d_sparse)),
    [], (0.0, 11.5); jac = true, sparse = true);
# BT.@btime DE.solve(prob_ode_brusselator_2d_mtk,Sundials.CVODE_BDF(linear_solver=:KLU),save_everystep=false); # compiles very slowly
nothing # hide
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
import LinearAlgebra
u0 = prob_ode_brusselator_2d_mtk.u0
p = prob_ode_brusselator_2d_mtk.p
const jaccache = prob_ode_brusselator_2d_mtk.f.jac(u0, p, 0.0)
const W = LinearAlgebra.I - 1.0 * jaccache

prectmp = IncompleteLU.ilu(W, τ = 50.0)
const preccache = Ref(prectmp)

function psetupilu(p, t, u, du, jok, jcurPtr, gamma)
    if jok
        prob_ode_brusselator_2d_mtk.f.jac(jaccache, u, p, t)
        jcurPtr[] = true

        # W = I - gamma*J
        @. W = -gamma * jaccache
        idxs = LinearAlgebra.diagind(W)
        @. @view(W[idxs]) = @view(W[idxs]) + 1

        # Build preconditioner on W
        preccache[] = IncompleteLU.ilu(W, τ = 5.0)
    end
end
```

Then the preconditioner action is to simply use the `ldiv!` of the generated
preconditioner:

```julia
function precilu(z, r, p, t, y, fy, gamma, delta, lr)
    ldiv!(z, preccache[], r)
end
```

We then simply pass these functions to the Sundials solver, with a choice of
`prec_side=1` to indicate that it is a left-preconditioner:

```julia
BT.@btime DE.solve(prob_ode_brusselator_2d_sparse,
    Sundials.CVODE_BDF(; linear_solver = :GMRES, prec = precilu, psetup = psetupilu,
        prec_side = 1); save_everystep = false);
```

And similarly for algebraic multigrid:

```julia
prectmp2 = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(W,
    presmoother = AlgebraicMultigrid.Jacobi(rand(size(W,
        1))),
    postsmoother = AlgebraicMultigrid.Jacobi(rand(size(W,
        1)))))
const preccache2 = Ref(prectmp2)
function psetupamg(p, t, u, du, jok, jcurPtr, gamma)
    if jok
        prob_ode_brusselator_2d_mtk.f.jac(jaccache, u, p, t)
        jcurPtr[] = true

        # W = I - gamma*J
        @. W = -gamma * jaccache
        idxs = LinearAlgebra.diagind(W)
        @. @view(W[idxs]) = @view(W[idxs]) + 1

        # Build preconditioner on W
        preccache2[] = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(W,
            presmoother = AlgebraicMultigrid.Jacobi(rand(size(W,
                1))),
            postsmoother = AlgebraicMultigrid.Jacobi(rand(size(W,
                1)))))
    end
end

function precamg(z, r, p, t, y, fy, gamma, delta, lr)
    ldiv!(z, preccache2[], r)
end

BT.@btime DE.solve(prob_ode_brusselator_2d_sparse,
    Sundials.CVODE_BDF(; linear_solver = :GMRES, prec = precamg, psetup = psetupamg,
        prec_side = 1); save_everystep = false);
```
