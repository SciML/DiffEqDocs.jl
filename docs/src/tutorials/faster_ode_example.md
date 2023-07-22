# [Code Optimization for Differential Equations](@id speed)

!!! note
    
    See [this FAQ](@ref faq_performance)
    for information on common pitfalls and how to improve performance.

## Code Optimization in Julia

Before starting this tutorial, we recommend the reader to check out one of the
many tutorials for optimization Julia code. The following is an incomplete
list:

  - [The Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
  - [MIT 18.337 Course Notes on Optimizing Serial Code](https://mitmath.github.io/18337/lecture2/optimizing)
  - [What scientists must know about hardware to write fast code](https://biojulia.dev/post/hardware/)

User-side optimizations are important because, for sufficiently difficult problems,
most time will be spent inside your `f` function, the function you are
trying to solve. “Efficient” integrators are those that reduce the required
number of `f` calls to hit the error tolerance. The main ideas for optimizing
your DiffEq code, or any Julia function, are the following:

  - Make it non-allocating
  - Use StaticArrays for small arrays
  - Use broadcast fusion
  - Make it type-stable
  - Reduce redundant calculations
  - Make use of BLAS calls
  - Optimize algorithm choice

We'll discuss these strategies in the context of differential equations.
Let's start with small systems.

## Example Accelerating a Non-Stiff Equation: The Lorenz Equation

Let's take the classic Lorenz system. Let's start by naively writing the
system in its out-of-place form:

```@example faster_ode
function lorenz(u, p, t)
    dx = 10.0 * (u[2] - u[1])
    dy = u[1] * (28.0 - u[3]) - u[2]
    dz = u[1] * u[2] - (8 / 3) * u[3]
    [dx, dy, dz]
end
```

Here, `lorenz` returns an object, `[dx,dy,dz]`, which is created within the body
of `lorenz`.

This is a common code pattern from high-level languages like MATLAB, SciPy, or
R's deSolve. However, the issue with this form is that it allocates a vector,
`[dx,dy,dz]`, at each step. Let's benchmark the solution process with this
choice of function:

```@example faster_ode
using DifferentialEquations, BenchmarkTools
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
@btime solve(prob, Tsit5());
nothing # hide
```

The `BenchmarkTools.jl` package's `@benchmark` runs the code multiple times to
get an accurate measurement. The minimum time is the time it takes when your
OS and other background processes aren't getting in the way. Notice that in
this case it takes about 5ms to solve and allocates around 11.11 MiB. However,
if we were to use this inside of a real user code, we'd see a lot of time spent
doing garbage collection (GC) to clean up all the arrays we made. Even if we
turn off saving, we have these allocations.

```@example faster_ode
@btime solve(prob, Tsit5(), save_everystep = false);
nothing # hide
```

The problem, of course, is that arrays are created every time our derivative
function is called. This function is called multiple times per step and is thus
the main source of memory usage. To fix this, we can use the in-place form to
***make our code non-allocating***:

```@example faster_ode
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    nothing
end
```

Here, instead of creating an array each time, we utilized the cache array `du`.
When the in-place form is used, DifferentialEquations.jl takes a different
internal route that minimizes the internal allocations as well.

!!! note
    
    Notice that nothing is returned. When in in-place form, the ODE solver ignores
    the return. Instead, make sure that the original `du` array is mutated instead
    of constructing a new array

When we benchmark this function, we will see quite a difference.

```@example faster_ode
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)
@btime solve(prob, Tsit5());
nothing # hide
```

```@example faster_ode
@btime solve(prob, Tsit5(), save_everystep = false);
nothing # hide
```

There is a 16x time difference just from that change! Notice there are still
some allocations and this is due to the construction of the integration cache.
But this doesn't scale with the problem size:

```@example faster_ode
tspan = (0.0, 500.0) # 5x longer than before
prob = ODEProblem(lorenz!, u0, tspan)
@btime solve(prob, Tsit5(), save_everystep = false);
nothing # hide
```

Since that's all setup allocations, the user-side optimization is complete.

### Further Optimizations of Small Non-Stiff ODEs with StaticArrays

Allocations are only expensive if they are “heap allocations”. For a more
in-depth definition of heap allocations,
[there are many sources online](http://net-informations.com/faq/net/stack-heap.htm).
But a good working definition is that heap allocations are variable-sized slabs
of memory which have to be pointed to, and this pointer indirection costs time.
Additionally, the heap has to be managed, and the garbage controllers has to
actively keep track of what's on the heap.

However, there's an alternative to heap allocations, known as stack allocations.
The stack is statically-sized (known at compile time) and thus its accesses are
quick. Additionally, the exact block of memory is known in advance by the
compiler, and thus re-using the memory is cheap. This means that allocating on
the stack has essentially no cost!

Arrays have to be heap allocated because their size (and thus the amount of
memory they take up) is determined at runtime. But there are structures in
Julia which are stack-allocated. `struct`s for example are stack-allocated
“value-type”s. `Tuple`s are a stack-allocated collection. The most useful data
structure for DiffEq though is the `StaticArray` from the package
[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). These arrays
have their length determined at compile-time. They are created using macros
attached to normal array expressions, for example:

```@example faster_ode
using StaticArrays
A = SA[2.0, 3.0, 5.0]
typeof(A) # SVector{3, Float64} (alias for SArray{Tuple{3}, Float64, 1, 3})
```

Notice that the `3` after `SVector` gives the size of the `SVector`. It cannot
be changed. Additionally, `SVector`s are immutable, so we have to create a new
`SVector` to change values. But remember, we don't have to worry about
allocations because this data structure is stack-allocated. `SArray`s have
numerous extra optimizations as well: they have fast matrix multiplication,
fast QR factorizations, etc. which directly make use of the information about
the size of the array. Thus, when possible, they should be used.

Unfortunately, static arrays can only be used for sufficiently small arrays.
After a certain size, they are forced to heap allocate after some instructions
and their compile time balloons. Thus, static arrays shouldn't be used if your
system has more than ~20 variables. Additionally, only the native Julia
algorithms can fully utilize static arrays.

Let's ***optimize `lorenz` using static arrays***. Note that in this case, we
want to use the out-of-place allocating form, but this time we want to output
a static array:

```@example faster_ode
function lorenz_static(u, p, t)
    dx = 10.0 * (u[2] - u[1])
    dy = u[1] * (28.0 - u[3]) - u[2]
    dz = u[1] * u[2] - (8 / 3) * u[3]
    SA[dx, dy, dz]
end
```

To make the solver internally use static arrays, we simply give it a static
array as the initial condition:

```@example faster_ode
u0 = SA[1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz_static, u0, tspan)
@btime solve(prob, Tsit5());
nothing # hide
```

```@example faster_ode
@btime solve(prob, Tsit5(), save_everystep = false);
nothing # hide
```

And that's pretty much all there is to it. With static arrays, you don't have to
worry about allocating, so use operations like `*` and don't worry about fusing
operations (discussed in the next section). Do “the vectorized code” of
R/MATLAB/Python and your code in this case will be fast, or directly use the
numbers/values.

## Example Accelerating a Stiff Equation: the Robertson Equation

For these next examples, let's solve the Robertson equations (also known as
ROBER):

```math
\begin{aligned}
\frac{dy_1}{dt} &= -0.04y₁ + 10^4 y_2 y_3 \\
\frac{dy_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
\frac{dy_3}{dt} &= 3*10^7 y_{2}^2 \\
\end{aligned}
```

Given that these equations are stiff, non-stiff ODE solvers like `Tsit5` or
`Vern9` will fail to solve these equations. The automatic algorithm will detect
this and automatically switch to something more robust to handle these issues.
For example:

```@example faster_ode2
using DifferentialEquations
using Plots
function rober!(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    nothing
end
prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1e5), [0.04, 3e7, 1e4])
sol = solve(prob)
plot(sol, tspan = (1e-2, 1e5), xscale = :log10)
```

```@example faster_ode2
using BenchmarkTools
@btime solve(prob);
nothing # hide
```

### Choosing a Good Solver

Choosing a good solver is required for getting top-notch speed. General
recommendations can be found on the solver page (for example, the
[ODE Solver Recommendations](@ref ode_solve)).
The current recommendations can be simplified to a Rosenbrock method
(`Rosenbrock23` or `Rodas5`) for smaller (<50 ODEs) problems, ESDIRK methods
for slightly larger (`TRBDF2` or `KenCarp4` for <2000 ODEs), and `QNDF` for even
larger problems. `lsoda` from [LSODA.jl](https://github.com/rveltz/LSODA.jl) is
sometimes worth a try for the medium-sized category.

More details on the solver to choose can be found by benchmarking. See the
[SciMLBenchmarks](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/) to
compare many solvers on many problems.

From this, we try the recommendation of `Rosenbrock23()` for stiff ODEs at
default tolerances:

```@example faster_ode2
@btime solve(prob, Rosenbrock23());
nothing # hide
```

### Declaring Jacobian Functions

In order to reduce the Jacobian construction cost, one can describe a Jacobian
function by using the `jac` argument for the `ODEFunction`. First we have to
derive the Jacobian ``\frac{df_i}{du_j}`` which is `J[i,j]`. From this, we get:

```@example faster_ode2
function rober_jac!(J, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    J[1, 1] = k₁ * -1
    J[2, 1] = k₁
    J[3, 1] = 0
    J[1, 2] = y₃ * k₃
    J[2, 2] = y₂ * k₂ * -2 + y₃ * k₃ * -1
    J[3, 2] = y₂ * 2 * k₂
    J[1, 3] = k₃ * y₂
    J[2, 3] = k₃ * y₂ * -1
    J[3, 3] = 0
    nothing
end
f! = ODEFunction(rober!, jac = rober_jac!)
prob_jac = ODEProblem(f!, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
```

```@example faster_ode2
@btime solve(prob_jac, Rosenbrock23());
nothing # hide
```

### Automatic Derivation of Jacobian Functions

But that was hard! If you want to take the symbolic Jacobian of numerical
code, we can make use of [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)
to symbolic-ify the numerical code and do the symbolic calculation and return
the Julia code for this.

```@example faster_ode2
using ModelingToolkit
de = modelingtoolkitize(prob)
```

We can tell it to compute the Jacobian if we want to see the code:

```@example faster_ode2
ModelingToolkit.generate_jacobian(de)[2] # Second is in-place
```

Now let's use that to give the analytical solution Jacobian:

```@example faster_ode2
prob_jac2 = ODEProblem(de, [], (0.0, 1e5), jac = true)
```

```@example faster_ode2
@btime solve(prob_jac2);
nothing # hide
```

See the [ModelingToolkit.jl documentation](https://docs.sciml.ai/ModelingToolkit/stable/) for more
details.

### Accelerating Small ODE Solves with Static Arrays

If the ODE is sufficiently small (<20 ODEs or so), using [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)
for the state variables can greatly enhance the performance. This is done by
making `u0` a `StaticArray` and writing an out-of-place non-mutating dispatch
for static arrays, for the ROBER problem, this looks like:

```@example faster_ode2
using StaticArrays
function rober_static(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du3 = k₂ * y₂^2
    SA[du1, du2, du3]
end
prob = ODEProblem(rober_static, SA[1.0, 0.0, 0.0], (0.0, 1e5), SA[0.04, 3e7, 1e4])
sol = solve(prob, Rosenbrock23())
```

If we benchmark this, we see a really fast solution with really low allocation
counts:

```@example faster_ode2
@btime sol = solve(prob, Rosenbrock23());
nothing # hide
```

This version is thus very amenable to multithreading and other forms of parallelism.

## Example Accelerating Linear Algebra PDE Semi-Discretization

In this tutorial, we will optimize the right-hand side definition of a PDE
semi-discretization.

!!! note
    
    We highly recommend looking at the [Solving Large Stiff Equations](@ref stiff)
    tutorial for details on customizing DifferentialEquations.jl for more
    efficient large-scale stiff ODE solving. This section will only focus on the
    user-side code.

Let's optimize the solution of a Reaction-Diffusion PDE's discretization.
In its discretized form, this is the ODE:

```math
\begin{align}
du &= D_1 (A_y u + u A_x) + \frac{au^2}{v} + \bar{u} - \alpha u\\
dv &= D_2 (A_y v + v A_x) + a u^2 + \beta v
\end{align}
```

where ``u``, ``v``, and ``A`` are matrices. Here, we will use the simplified
version where ``A`` is the tridiagonal stencil ``[1,-2,1]``, i.e. it's the 2D
discretization of the Laplacian. The native code would be something along the
lines of:

```@example faster_ode3
using DifferentialEquations, LinearAlgebra, BenchmarkTools
# Generate the constants
p = (1.0, 1.0, 1.0, 10.0, 0.001, 100.0) # a,α,ubar,β,D1,D2
N = 100
Ax = Array(Tridiagonal([1.0 for i in 1:(N - 1)], [-2.0 for i in 1:N],
    [1.0 for i in 1:(N - 1)]))
Ay = copy(Ax)
Ax[2, 1] = 2.0
Ax[end - 1, end] = 2.0
Ay[1, 2] = 2.0
Ay[end, end - 1] = 2.0

function basic_version!(dr, r, p, t)
    a, α, ubar, β, D1, D2 = p
    u = r[:, :, 1]
    v = r[:, :, 2]
    Du = D1 * (Ay * u + u * Ax)
    Dv = D2 * (Ay * v + v * Ax)
    dr[:, :, 1] = Du .+ a .* u .* u ./ v .+ ubar .- α * u
    dr[:, :, 2] = Dv .+ a .* u .* u .- β * v
end

a, α, ubar, β, D1, D2 = p
uss = (ubar + β) / α
vss = (a / β) * uss^2
r0 = zeros(100, 100, 2)
r0[:, :, 1] .= uss .+ 0.1 .* rand.()
r0[:, :, 2] .= vss

prob = ODEProblem(basic_version!, r0, (0.0, 0.1), p)
```

In this version, we have encoded our initial condition to be a 3-dimensional
array, with `u[:,:,1]` being the `A` part and `u[:,:,2]` being the `B` part.

```@example faster_ode3
@btime solve(prob, Tsit5());
nothing # hide
```

While this version isn't very efficient,

#### We recommend writing the “high-level” code first, and iteratively optimizing it!

The first thing that we can do is get rid of the slicing allocations. The
operation `r[:,:,1]` creates a temporary array instead of a “view”, i.e. a
pointer to the already existing memory. To make it a view, add `@view`. Note
that we have to be careful with views because they point to the same memory,
and thus changing a view changes the original values:

```@example faster_ode3
A = rand(4)
@show A
B = @view A[1:3]
B[2] = 2
@show A
```

Notice that changing `B` changed `A`. This is something to be careful of, but
at the same time we want to use this since we want to modify the output `dr`.
Additionally, the last statement is a purely element-wise operation, and thus
we can make use of broadcast fusion there. Let's rewrite `basic_version!` to
***avoid slicing allocations*** and to ***use broadcast fusion***:

```@example faster_ode3
function gm2!(dr, r, p, t)
    a, α, ubar, β, D1, D2 = p
    u = @view r[:, :, 1]
    v = @view r[:, :, 2]
    du = @view dr[:, :, 1]
    dv = @view dr[:, :, 2]
    Du = D1 * (Ay * u + u * Ax)
    Dv = D2 * (Ay * v + v * Ax)
    @. du = Du + a .* u .* u ./ v + ubar - α * u
    @. dv = Dv + a .* u .* u - β * v
end
prob = ODEProblem(gm2!, r0, (0.0, 0.1), p)
@btime solve(prob, Tsit5());
nothing # hide
```

Now, most of the allocations are taking place in `Du = D1*(Ay*u + u*Ax)` since
those operations are vectorized and not mutating. We should instead replace the
matrix multiplications with `mul!`. When doing so, we will need to have cache
variables to write into. This looks like:

```@example faster_ode3
Ayu = zeros(N, N)
uAx = zeros(N, N)
Du = zeros(N, N)
Ayv = zeros(N, N)
vAx = zeros(N, N)
Dv = zeros(N, N)
function gm3!(dr, r, p, t)
    a, α, ubar, β, D1, D2 = p
    u = @view r[:, :, 1]
    v = @view r[:, :, 2]
    du = @view dr[:, :, 1]
    dv = @view dr[:, :, 2]
    mul!(Ayu, Ay, u)
    mul!(uAx, u, Ax)
    mul!(Ayv, Ay, v)
    mul!(vAx, v, Ax)
    @. Du = D1 * (Ayu + uAx)
    @. Dv = D2 * (Ayv + vAx)
    @. du = Du + a * u * u ./ v + ubar - α * u
    @. dv = Dv + a * u * u - β * v
end
prob = ODEProblem(gm3!, r0, (0.0, 0.1), p)
@btime solve(prob, Tsit5());
nothing # hide
```

But our temporary variables are global variables. We need to either declare the caches as `const` or localize them. We can localize them by adding them to the parameters, `p`. It's easier for the compiler to reason about local variables than global variables. ***Localizing variables helps to ensure type stability***.

```@example faster_ode3
p = (1.0, 1.0, 1.0, 10.0, 0.001, 100.0, Ayu, uAx, Du, Ayv, vAx, Dv) # a,α,ubar,β,D1,D2
function gm4!(dr, r, p, t)
    a, α, ubar, β, D1, D2, Ayu, uAx, Du, Ayv, vAx, Dv = p
    u = @view r[:, :, 1]
    v = @view r[:, :, 2]
    du = @view dr[:, :, 1]
    dv = @view dr[:, :, 2]
    mul!(Ayu, Ay, u)
    mul!(uAx, u, Ax)
    mul!(Ayv, Ay, v)
    mul!(vAx, v, Ax)
    @. Du = D1 * (Ayu + uAx)
    @. Dv = D2 * (Ayv + vAx)
    @. du = Du + a * u * u ./ v + ubar - α * u
    @. dv = Dv + a * u * u - β * v
end
prob = ODEProblem(gm4!, r0, (0.0, 0.1), p)
@btime solve(prob, Tsit5());
nothing # hide
```

We could then use the BLAS `gemmv` to optimize the matrix multiplications some more, but instead let's devectorize the stencil.

```@example faster_ode3
p = (1.0, 1.0, 1.0, 10.0, 0.001, 100.0, N)
function fast_gm!(du, u, p, t)
    a, α, ubar, β, D1, D2, N = p

    @inbounds for j in 2:(N - 1), i in 2:(N - 1)
        du[i, j, 1] = D1 *
                      (u[i - 1, j, 1] + u[i + 1, j, 1] + u[i, j + 1, 1] + u[i, j - 1, 1] -
                       4u[i, j, 1]) +
                      a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
    end

    @inbounds for j in 2:(N - 1), i in 2:(N - 1)
        du[i, j, 2] = D2 *
                      (u[i - 1, j, 2] + u[i + 1, j, 2] + u[i, j + 1, 2] + u[i, j - 1, 2] -
                       4u[i, j, 2]) +
                      a * u[i, j, 1]^2 - β * u[i, j, 2]
    end

    @inbounds for j in 2:(N - 1)
        i = 1
        du[1, j, 1] = D1 *
                      (2u[i + 1, j, 1] + u[i, j + 1, 1] + u[i, j - 1, 1] - 4u[i, j, 1]) +
                      a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
    end
    @inbounds for j in 2:(N - 1)
        i = 1
        du[1, j, 2] = D2 *
                      (2u[i + 1, j, 2] + u[i, j + 1, 2] + u[i, j - 1, 2] - 4u[i, j, 2]) +
                      a * u[i, j, 1]^2 - β * u[i, j, 2]
    end
    @inbounds for j in 2:(N - 1)
        i = N
        du[end, j, 1] = D1 *
                        (2u[i - 1, j, 1] + u[i, j + 1, 1] + u[i, j - 1, 1] - 4u[i, j, 1]) +
                        a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
    end
    @inbounds for j in 2:(N - 1)
        i = N
        du[end, j, 2] = D2 *
                        (2u[i - 1, j, 2] + u[i, j + 1, 2] + u[i, j - 1, 2] - 4u[i, j, 2]) +
                        a * u[i, j, 1]^2 - β * u[i, j, 2]
    end

    @inbounds for i in 2:(N - 1)
        j = 1
        du[i, 1, 1] = D1 *
                      (u[i - 1, j, 1] + u[i + 1, j, 1] + 2u[i, j + 1, 1] - 4u[i, j, 1]) +
                      a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
    end
    @inbounds for i in 2:(N - 1)
        j = 1
        du[i, 1, 2] = D2 *
                      (u[i - 1, j, 2] + u[i + 1, j, 2] + 2u[i, j + 1, 2] - 4u[i, j, 2]) +
                      a * u[i, j, 1]^2 - β * u[i, j, 2]
    end
    @inbounds for i in 2:(N - 1)
        j = N
        du[i, end, 1] = D1 *
                        (u[i - 1, j, 1] + u[i + 1, j, 1] + 2u[i, j - 1, 1] - 4u[i, j, 1]) +
                        a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
    end
    @inbounds for i in 2:(N - 1)
        j = N
        du[i, end, 2] = D2 *
                        (u[i - 1, j, 2] + u[i + 1, j, 2] + 2u[i, j - 1, 2] - 4u[i, j, 2]) +
                        a * u[i, j, 1]^2 - β * u[i, j, 2]
    end

    @inbounds begin
        i = 1
        j = 1
        du[1, 1, 1] = D1 * (2u[i + 1, j, 1] + 2u[i, j + 1, 1] - 4u[i, j, 1]) +
                      a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
        du[1, 1, 2] = D2 * (2u[i + 1, j, 2] + 2u[i, j + 1, 2] - 4u[i, j, 2]) +
                      a * u[i, j, 1]^2 - β * u[i, j, 2]

        i = 1
        j = N
        du[1, N, 1] = D1 * (2u[i + 1, j, 1] + 2u[i, j - 1, 1] - 4u[i, j, 1]) +
                      a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
        du[1, N, 2] = D2 * (2u[i + 1, j, 2] + 2u[i, j - 1, 2] - 4u[i, j, 2]) +
                      a * u[i, j, 1]^2 - β * u[i, j, 2]

        i = N
        j = 1
        du[N, 1, 1] = D1 * (2u[i - 1, j, 1] + 2u[i, j + 1, 1] - 4u[i, j, 1]) +
                      a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
        du[N, 1, 2] = D2 * (2u[i - 1, j, 2] + 2u[i, j + 1, 2] - 4u[i, j, 2]) +
                      a * u[i, j, 1]^2 - β * u[i, j, 2]

        i = N
        j = N
        du[end, end, 1] = D1 * (2u[i - 1, j, 1] + 2u[i, j - 1, 1] - 4u[i, j, 1]) +
                          a * u[i, j, 1]^2 / u[i, j, 2] + ubar - α * u[i, j, 1]
        du[end, end, 2] = D2 * (2u[i - 1, j, 2] + 2u[i, j - 1, 2] - 4u[i, j, 2]) +
                          a * u[i, j, 1]^2 - β * u[i, j, 2]
    end
end
prob = ODEProblem(fast_gm!, r0, (0.0, 0.1), p)
@btime solve(prob, Tsit5());
nothing # hide
```

Notice that in this case fusing the loops and avoiding the linear operators is a
major improvement of about 10x! That's an order of magnitude faster than our
original MATLAB/SciPy/R vectorized style code!

Since this is tedious to do by hand, we note that
[ModelingToolkit.jl's symbolic code generation](https://docs.sciml.ai/ModelingToolkit/stable/) can
do this automatically from the basic version:

```@example faster_ode3
using ModelingToolkit
function basic_version!(dr, r, p, t)
    a, α, ubar, β, D1, D2 = p
    u = r[:, :, 1]
    v = r[:, :, 2]
    Du = D1 * (Ay * u + u * Ax)
    Dv = D2 * (Ay * v + v * Ax)
    dr[:, :, 1] = Du .+ a .* u .* u ./ v .+ ubar .- α * u
    dr[:, :, 2] = Dv .+ a .* u .* u .- β * v
end

a, α, ubar, β, D1, D2 = p
uss = (ubar + β) / α
vss = (a / β) * uss^2
r0 = zeros(100, 100, 2)
r0[:, :, 1] .= uss .+ 0.1 .* rand.()
r0[:, :, 2] .= vss

prob = ODEProblem(basic_version!, r0, (0.0, 0.1), p)
de = modelingtoolkitize(prob)

# Note jac=true,sparse=true makes it automatically build sparse Jacobian code
# as well!

fastprob = ODEProblem(de, [], (0.0, 0.1), jac = true, sparse = true)
```

Lastly, we can do other things like multithread the main loops.
[LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) provides
the `@turbo` macro for doing a lot of SIMD enhancements, and `@tturbo` is the
multithreaded version.

### Optimizing Algorithm Choices

The last thing to do is then ***optimize our algorithm choice***. We have been
using `Tsit5()` as our test algorithm, but in reality this problem is a stiff
PDE discretization and thus one recommendation is to use `CVODE_BDF()`. However,
instead of using the default dense Jacobian, we should make use of the sparse
Jacobian afforded by the problem. The Jacobian is the matrix $\frac{df_i}{dr_j}$,
where $r$ is read by the linear index (i.e. down columns). But since the $u$
variables depend on the $v$, the band size here is large, and thus this will
not do well with a Banded Jacobian solver. Instead, we utilize sparse Jacobian
algorithms. `CVODE_BDF` allows us to use a sparse Newton-Krylov solver by
setting `linear_solver = :GMRES`.

!!! note
    
    The [Solving Large Stiff Equations](@ref stiff) tutorial goes through these
    details. This is simply to give a taste of how much optimization opportunity
    is left on the table!

Let's see how our fast right-hand side scales as we increase the integration time.

```@example faster_ode3
prob = ODEProblem(fast_gm!, r0, (0.0, 10.0), p)
@btime solve(prob, Tsit5());
nothing # hide
```

```@example faster_ode3
using Sundials
@btime solve(prob, CVODE_BDF(linear_solver = :GMRES));
nothing # hide
```

```@example faster_ode3
prob = ODEProblem(fast_gm!, r0, (0.0, 100.0), p)
# Will go out of memory if we don't turn off `save_everystep`!
@btime solve(prob, Tsit5(), save_everystep = false);
nothing # hide
```

```@example faster_ode3
@btime solve(prob, CVODE_BDF(linear_solver = :GMRES), save_everystep = false);
nothing # hide
```

```@example faster_ode3
prob = ODEProblem(fast_gm!, r0, (0.0, 500.0), p)
@btime solve(prob, CVODE_BDF(linear_solver = :GMRES), save_everystep = false);
nothing # hide
```

Notice that we've eliminated almost all allocations, allowing the code to grow
without hitting garbage collection and slowing down.

Why is `CVODE_BDF` doing well? What's happening is that, because the problem is
stiff, the number of steps required by the explicit Runge-Kutta method grows
rapidly, whereas `CVODE_BDF` is taking large steps. Additionally, the `GMRES`
linear solver form is quite an efficient way to solve the implicit system in
this case. This is problem-dependent, and in many cases using a Krylov method
effectively requires a preconditioner, so you need to play around with testing
other algorithms and linear solvers to find out what works best with your
problem.

Now continue to the [Solving Large Stiff Equations](@ref stiff) tutorial for more
details on optimizing the algorithm choice for such codes.
