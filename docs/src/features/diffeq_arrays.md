# DiffEq-Specific Array Types

In many cases, a standard array may not be enough to fully hold the data for a
model. Many of the solvers in DifferentialEquations.jl (only the native Julia
methods) allow you to solve problems on `AbstractArray` types which allow you
to extend the meaning of an array. This page describes some of the `AbstractArray`
types which can be helpful for modeling differential equations problems.

## ArrayPartitions

ArrayPartitions in DiffEq are used for heterogeneous arrays. For example,
`PartitionedODEProblem` solvers use them internally to turn the separate parts
into a single array. You can construct an `ArrayPartition` using RecursiveArrayTools.jl:

```julia
using RecursiveArrayTools
A = ArrayPartition(x::AbstractArray...)
```

where is `x` a list of arrays. The resulting `A` will act like a single array, and its
broadcast will be type stable, allowing for it to be used inside of the native Julia
DiffEq solvers in an efficient way. This is a good way to generate an array which
has different units for different parts, or different amounts of precision.

### Usage

An `ArrayPartition` acts like a single array. `A[i]` indexes through the first
array, then the second, etc. all linearly. But `A.x` is where the arrays are stored.
Thus for

```julia
using RecursiveArrayTools
A = ArrayPartition(y,z)
```

We would have `A.x[1]==y` and `A.x[2]==z`. Broadcasting like `f.(A)` is efficient.

### Example: Dynamics Equations

In this example we will show using heterogeneous units in dynamics equations. Our
arrays will be:

```julia
using Unitful, RecursiveArrayTools, DiffEqBase, OrdinaryDiffEq
using LinearAlgebra

r0 = [1131.340, -2282.343, 6672.423]u"km"
v0 = [-5.64305, 4.30333, 2.42879]u"km/s"
Δt = 86400.0*365u"s"
μ = 398600.4418u"km^3/s^2"
rv0 = ArrayPartition(r0,v0)
```

Here, `r0` is the initial positions, and `v0` are the initial velocities. `rv0`
is the `ArrayPartition` initial condition. We now write our update function in
terms of the `ArrayPartition`:

```julia
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

```julia
prob = ODEProblem(f, rv0, (0.0u"s", Δt), μ)
sol = solve(prob, Vern8())
```

## MultiScaleArrays

The multi-scale modeling functionality is provided by MultiScaleArrays.jl. It
allows for designing a multi-scale model as an extension of an array, which in
turn can be directly used in the native Julia solvers of DifferentialEquations.jl.

For more information, please see [the MultiScaleArrays.jl README](https://github.com/JuliaDiffEq/MultiScaleArrays.jl).

## DEDataArrays

The `DEDataArray{T}` type allows one to add other "non-continuous" variables
to an array, which can be useful in many modeling situations involving lots of
events. To define an `DEDataArray`, make a type which subtypes `DEDataArray{T}`
with a field `x` for the "array of continuous variables" for which you would
like the differential equation to treat directly. The other fields are treated
as "discrete variables". For example:

```julia
mutable struct MyDataArray{T,N} <: DEDataArray{T,N}
    x::Array{T,1}
    a::T
    b::Symbol
end
```

In this example, our resultant array is a `SimType`, and its data which is presented
to the differential equation solver will be the array `x`. Any array which the
differential equation solver can use is allowed to be made as the field `x`, including
other `DEDataArray`s. Other than that, you can add whatever fields you please, and
let them be whatever type you please.

These extra fields are carried along in the differential equation solver that
the user can use in their `f` equation and modify via callbacks. For example,
inside of a an update function, it is safe to do:

```julia
function f(du,u,p,t)
  u.a = t
end
```

to update the discrete variables (unless the algorithm notes that it does not
step to the endpoint, in which case a callback must be used to update appropriately.)

Note that the aliases `DEDataVector` and `DEDataMatrix` cover the one and two
dimensional cases.

### Example: A Control Problem

In this example we will use a `DEDataArray` to solve a problem where control parameters
change at various timepoints. First we will build

```julia
mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end
```

as our `DEDataVector`. It has an extra field `f1` which we will use as our control
variable. Our ODE function will use this field as follows:

```julia
function f(du,u,p,t)
    du[1] = -0.5*u[1] + u.f1
    du[2] = -0.5*u[2]
end
```

Now we will setup our control mechanism. It will be a simple setup which uses
set timepoints at which we will change `f1`. At `t=5.0` we will want to increase
the value of `f1`, and at `t=8.0` we will want to decrease the value of `f1`. Using
the [`DiscreteCallback` interface](../callback_functions.html), we code these conditions
as follows:

```julia
const tstop1 = [5.]
const tstop2 = [8.]


function condition(u,t,integrator)
  t in tstop1
end

function condition2(u,t,integrator)
  t in tstop2
end
```

Now we have to apply an effect when these conditions are reached. When `condition`
is hit (at `t=5.0`), we will increase `f1` to 1.5. When `condition2` is reached,
we will decrease `f1` to `-1.5`. This is done via the functions:

```julia
function affect!(integrator)
  for c in full_cache(integrator)
    c.f1 = 1.5
  end
end

function affect2!(integrator)
  for c in full_cache(integrator)
    c.f1 = -1.5
  end
end
```

Notice that we have to loop through the `full_cache` array (provided by the integrator
interface) to ensure that all internal caches are also updated. With these functions
we can build our callbacks:

```julia
save_positions = (true,true)

cb = DiscreteCallback(condition, affect!, save_positions=save_positions)

save_positions = (false,true)

cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)

cbs = CallbackSet(cb,cb2)
```


Now we define our initial condition. We will start at `[10.0;10.0]` with `f1=0.0`.

```julia
u0 = SimType([10.0;10.0], 0.0)
prob = ODEProblem(f,u0,(0.0,10.0))
```

Lastly we solve the problem. Note that we must pass `tstop` values of `5.0` and
`8.0` to ensure the solver hits those timepoints exactly:

```julia
const tstop = [5.;8.]
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
```

![data_array_plot](../assets/data_array.png)

It's clear from the plot how the controls affected the outcome.

### Data Arrays vs ParameterizedFunctions

The reason for using a `DEDataArray` is because the solution will then save the
control parameters. For example, we can see what the control parameter was at
every timepoint by checking:

```julia
[sol[i].f1 for i in 1:length(sol)]
```

A similar solution can be achieved using a `ParameterizedFunction`.
We could have instead created our function as:

```julia
function f(du,u,p,t)
    du[1] = -0.5*u[1] + p
    du[2] = -0.5*u[2]
end
u0 = SimType([10.0;10.0], 0.0)
p = 0.0
prob = ODEProblem(f,u0,(0.0,10.0),p)
const tstop = [5.;8.]
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
```

where we now change the callbacks to changing the parameter:

```julia
function affect!(integrator)
  integrator.p = 1.5
end

function affect2!(integrator)
  integrator.p = -1.5
end
```

This will also solve the equation and get a similar result. It will also be slightly
faster in some cases. However, if the equation is solved in this manner, there will
be no record of what the parameter was at each timepoint. That is the tradeoff
between `DEDataArray`s and `ParameterizedFunction`s.
