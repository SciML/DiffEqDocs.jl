# [Code Optimization for Differential Equations](@ref speed)

!!! note

    See [this FAQ](@ref faq_performance)
    for information on common pitfalls and how to improve performance. For a similar
    tutorial, see [Optimizing DiffEq Code](http://tutorials.juliadiffeq.org/html/introduction/03-optimizing_diffeq_code.html).

## Code Optimization in Julia

Before starting this tutorial, we recommend the reader to check out one of the
many tutorials for optimization Julia code. The following is an incomplete
list:

- [The Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
- [MIT 18.337 Course Notes on Optimizing Serial Code](https://mitmath.github.io/18337/lecture2/optimizing)
- [What scientists must know about hardware to write fast code](https://biojulia.net/post/hardware/)

User-side optimizations are important because, for sufficiently difficult problems, most of the time will be spent inside of your `f` function, the function you are trying to solve. "Efficient" integrators are those that reduce the required number of `f` calls to hit the error tolerance. The main ideas for optimizing your DiffEq code, or any Julia function, are the following:

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

```julia
function lorenz(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 [dx,dy,dz]
end
```

Here, `lorenz` returns an object, `[dx,dy,dz]`, which is created within the body
of `lorenz`.

This is a common code pattern from high-level languages like MATLAB, SciPy, or
R's deSolve. However, the issue with this form is that it allocates a vector,
`[dx,dy,dz]`, at each step. Let's benchmark the solution process with this
choice of function:

```julia
using DifferentialEquations, BenchmarkTools
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
@benchmark solve(prob,Tsit5())
```

```julia
BenchmarkTools.Trial: 1350 samples with 1 evaluation.
 Range (min … max):  2.296 ms … 13.764 ms  ┊ GC (min … max):  0.00% … 67.48%
 Time  (median):     2.561 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   3.699 ms ±  2.223 ms  ┊ GC (mean ± σ):  14.83% ± 17.79%

  █▆▄       ▁ ▃▄▃▄▂
  ████▆▆█▇▇▇█▇█████▆▄▄▄▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▄▆▇▆▇▇▆▇▇▆█▇▇▆ █
  2.3 ms       Histogram: log(frequency) by time     11.3 ms <

 Memory estimate: 7.82 MiB, allocs estimate: 101102.
```

The `BenchmarkTools.jl` package's `@benchmark` runs the code multiple times to
get an accurate measurement. The minimum time is the time it takes when your
OS and other background processes aren't getting in the way. Notice that in
this case it takes about 5ms to solve and allocates around 11.11 MiB. However,
if we were to use this inside of a real user code we'd see a lot of time spent
doing garbage collection (GC) to clean up all of the arrays we made. Even if we
turn off saving we have these allocations.

```julia
@benchmark solve(prob,Tsit5(),save_everystep=false)

BenchmarkTools.Trial: 1490 samples with 1 evaluation.
 Range (min … max):  2.010 ms … 14.612 ms  ┊ GC (min … max):  0.00% … 65.56%
 Time  (median):     2.313 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   3.350 ms ±  2.095 ms  ┊ GC (mean ± σ):  14.55% ± 17.48%

  █▇▅▁       ▃▅▅▄▃
  █████▇▇▇▆▇███████▅▁▄▁▁▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▄▄▅▇▇▆▆█▇▇▆▆▆▆ █
  2.01 ms      Histogram: log(frequency) by time       11 ms <

 Memory estimate: 6.83 MiB, allocs estimate: 89529.
```

The problem of course is that arrays are created every time our derivative
function is called. This function is called multiple times per step and is thus
the main source of memory usage. To fix this, we can use the in-place form to
***make our code non-allocating***:

```julia
function lorenz!(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
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

```julia
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5())

BenchmarkTools.Trial: 8180 samples with 1 evaluation.
 Range (min … max):  415.800 μs …  12.112 ms  ┊ GC (min … max):  0.00% … 93.21%
 Time  (median):     463.700 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   605.847 μs ± 695.892 μs  ┊ GC (mean ± σ):  11.02% ±  9.07%

  ▄█▇▅▃▂▂▂▂▁▁▄▅▅▃▂▂▂▁▂▁▁                                        ▂
  ███████████████████████▇▇▇▇▇▆▆▆▅▅▄▄▅▃▅▄▅▄▅▄▁▅▁▃▄▁▃▁▁▁▁▃▄▁▁▄▁▃ █
  416 μs        Histogram: log(frequency) by time       1.68 ms <

 Memory estimate: 1016.36 KiB, allocs estimate: 11641.
```

```julia
@benchmark solve(prob,Tsit5(),save_everystep=false)

BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  197.900 μs … 315.800 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     206.700 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   207.688 μs ±   6.241 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

          ▃▅▅▇█▅▃       ▁▁▂▁▃▂▁
  ▁▁▁▂▃▃▅▇█████████▇▇▇▆█████████▆▆▅▅▃▂▃▂▂▂▂▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  198 μs           Histogram: frequency by time          227 μs <

 Memory estimate: 4.94 KiB, allocs estimate: 41.
```

There is a 16x time difference just from that change! Notice there are still
some allocations and this is due to the construction of the integration cache.
But this doesn't scale with the problem size:

```julia
tspan = (0.0,500.0) # 5x longer than before
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5(),save_everystep=false)

BenchmarkTools.Trial: 4755 samples with 1 evaluation.
 Range (min … max):  1.014 ms …  1.485 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     1.042 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.048 ms ± 31.281 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▁▆▆▅▇█▇▆▆▂▁ ▁            ▁
  ▂▄█████████████▇▇▇▇▇▇█▇▇█▇███▇█▇▇▅▆▄▄▄▃▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁ ▄
  1.01 ms        Histogram: frequency by time        1.12 ms <

 Memory estimate: 4.94 KiB, allocs estimate: 41.
```

Since that's all setup allocations the user-side optimization is complete.

### Further Optimizations of Small Non-Stiff ODEs with StaticArrays

Allocations are only expensive if they are "heap allocations". For a more
in-depth definition of heap allocations,
[there are a lot of sources online](http://net-informations.com/faq/net/stack-heap.htm).
But a good working definition is that heap allocations are variable-sized slabs
of memory which have to be pointed to, and this pointer indirection costs time.
Additionally, the heap has to be managed and the garbage controllers has to
actively keep track of what's on the heap.

However, there's an alternative to heap allocations, known as stack allocations.
The stack is statically-sized (known at compile time) and thus its accesses are
quick. Additionally, the exact block of memory is known in advance by the
compiler, and thus re-using the memory is cheap. This means that allocating on
the stack has essentially no cost!

Arrays have to be heap allocated because their size (and thus the amount of
memory they take up) is determined at runtime. But there are structures in
Julia which are stack-allocated. `struct`s for example are stack-allocated
"value-type"s. `Tuple`s are a stack-allocated collection. The most useful data
structure for DiffEq though is the `StaticArray` from the package
[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). These arrays
have their length determined at compile-time. They are created using macros
attached to normal array expressions, for example:

```julia
using StaticArrays
A = SA[2.0,3.0,5.0]
typeof(A) # SVector{3, Float64} (alias for SArray{Tuple{3}, Float64, 1, 3})
```

Notice that the `3` after `SVector` gives the size of the `SVector`. It cannot
be changed. Additionally, `SVector`s are immutable, so we have to create a new
`SVector` to change values. But remember, we don't have to worry about
allocations because this data structure is stack-allocated. `SArray`s have a
lot of extra optimizations as well: they have fast matrix multiplication,
fast QR factorizations, etc. which directly make use of the information about
the size of the array. Thus, when possible they should be used.

Unfortunately static arrays can only be used for sufficiently small arrays.
After a certain size, they are forced to heap allocate after some instructions
and their compile time balloons. Thus static arrays shouldn't be used if your
system has more than ~20 variables. Additionally, only the native Julia
algorithms can fully utilize static arrays.

Let's ***optimize `lorenz` using static arrays***. Note that in this case, we
want to use the out-of-place allocating form, but this time we want to output
a static array:

```julia
function lorenz_static(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 SA[dx,dy,dz]
end
```

To make the solver internally use static arrays, we simply give it a static array as the initial condition:

```julia
u0 = SA[1.0,0.0,0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz_static,u0,tspan)
@benchmark solve(prob,Tsit5())

BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  196.600 μs …   6.310 ms  ┊ GC (min … max): 0.00% … 95.70%
 Time  (median):     220.900 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   265.006 μs ± 302.623 μs  ┊ GC (mean ± σ):  6.91% ±  5.82%

    ▅█▄▄▄
  ▁▄██████▅▄▃▂▂▁▁▂▁▁▁▁▁▁▁▂▃▄▅▅▆▆▆▅▅▅▄▃▃▃▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  197 μs           Histogram: frequency by time          374 μs <

 Memory estimate: 394.50 KiB, allocs estimate: 1319.
```

```julia
@benchmark solve(prob,Tsit5(),save_everystep=false)

BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  144.100 μs … 242.600 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     151.000 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   151.875 μs ±   7.502 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

      █▇▇▂   ▁▅▄▂▂
  ▂▂▃█████▇▇▇██████▅▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▂ ▃
  144 μs           Histogram: frequency by time          185 μs <

 Memory estimate: 3.67 KiB, allocs estimate: 22.
```

And that's pretty much all there is to it. With static arrays you don't have to
worry about allocating, so use operations like `*` and don't worry about fusing
operations (discussed in the next section). Do "the vectorized code" of
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

Now let's improve that.

### Declaring Jacobian Functions

In order to reduce the Jacobian construction cost, one can describe a Jacobian
function by using the `jac` argument for the `ODEFunction`. First we have to
derive the Jacobian ``\frac{df_i}{du_j}`` which is `J[i,j]`. From this we get:

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
to symbolic-ify the numerical code and do the symbolic calculation and return
the Julia code for this.

```julia
using ModelingToolkit
de = modelingtoolkitize(prob)
```

We can tell it to compute the Jacobian if we want to see the code:

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
prob_jac2 = ODEProblem(de,[],(0.0,1e5),jac=true)
```

```julia-repl
julia> @btime solve(prob_jac2)
122.600 μs (1425 allocations: 99.34 KiB)
```

See the [ModelingToolkit.jl documentation](https://mtk.sciml.ai/dev/) for more
details.

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
