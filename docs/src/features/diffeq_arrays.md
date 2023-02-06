# [DiffEq-Specific Array Types](@id diffeq_arrays)

In many cases, a standard array may not be enough to fully hold the data for a
model. Many of the solvers in DifferentialEquations.jl (only the native Julia
methods) allow you to solve problems on `AbstractArray` types, which allow you
to extend the meaning of an array. This page describes some `AbstractArray`
types which can be helpful for modeling differential equations problems.

## ArrayPartitions

ArrayPartitions in DiffEq are used for heterogeneous arrays. For example,
`DynamicalODEProblem` solvers use them internally to turn the separate parts
into a single array. You can construct an `ArrayPartition` using RecursiveArrayTools.jl:

```julia
using RecursiveArrayTools
A = ArrayPartition(x::AbstractArray...)
```

where `x` is an array of arrays. Then, `A` will act like a single array, and its
broadcast will be type stable, allowing for it to be efficiently used inside
the native Julia DiffEq solvers. This is a good way to generate an array which
has different units for different parts, or different amounts of precision.

### Usage

An `ArrayPartition` acts like a single array. `A[i]` indexes through the first
array, then the second, etc. all linearly. But `A.x` is where the arrays are stored.
Thus for

```julia
using RecursiveArrayTools
A = ArrayPartition(y, z)
```

We would have `A.x[1]==y` and `A.x[2]==z`. Broadcasting like `f.(A)` is efficient.

### Example: Dynamics Equations

In this example, we will show using heterogeneous units in dynamics equations. Our
arrays will be:

```@example diffeq_arrays
using Unitful, RecursiveArrayTools, DifferentialEquations
using LinearAlgebra

r0 = [1131.340, -2282.343, 6672.423]u"km"
v0 = [-5.64305, 4.30333, 2.42879]u"km/s"
Δt = 86400.0 * 365u"s"
μ = 398600.4418u"km^3/s^2"
rv0 = ArrayPartition(r0, v0)
```

Here, `r0` is the initial positions, and `v0` are the initial velocities. `rv0`
is the `ArrayPartition` initial condition. We now write our update function in
terms of the `ArrayPartition`:

```@example diffeq_arrays
function f(dy, y, μ, t)
    r = norm(y.x[1])
    dy.x[1] .= y.x[2]
    dy.x[2] .= -μ .* y.x[1] / r^3
end
```

Notice that `y.x[1]` is the `r` part of `y`, and `y.x[2]` is the `v` part of `y`.
Using this kind of indexing is type stable, even though the array itself is
heterogeneous. Note that one can also use things like `2y` or `y.+x` and the
broadcasting will be efficient.

Now to solve our equations, we do the same thing as always in DiffEq:

```@example diffeq_arrays
prob = ODEProblem(f, rv0, (0.0u"s", Δt), μ)
sol = solve(prob, Vern8())
```

## MultiScaleArrays

The multi-scale modeling functionality is provided by MultiScaleArrays.jl. It
allows for designing a multi-scale model as an extension of an array, which in
turn can be directly used in the native Julia solvers of DifferentialEquations.jl.

For more information, please see [the MultiScaleArrays.jl README](https://github.com/SciML/MultiScaleArrays.jl).
