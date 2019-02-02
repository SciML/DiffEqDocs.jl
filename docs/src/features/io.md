# I/O: Saving and Loading Solution Data

The ability to save and load solutions is important for handling large datasets
and analyzing the results over multiple Julia sessions. This page explains the
existing functionality for doing so.

## Tabular Data: IterableTables

An interface to [IterableTables.jl](https://github.com/davidanthoff/IterableTables.jl)
is provided. This IterableTables link allows you to use a solution
type as the data source to convert to other tabular data formats. For example,
let's solve a 4x2 system of ODEs:

```julia
f_2dlinear = (du,u,p,t) -> du.=1.01u
prob = ODEProblem(f_2dlinear,rand(2,2),(0.0,1.0))
sol1 =solve(prob,Euler();dt=1//2^(4))
```

then we can convert this to a dataframe using `DataFrame`:

```julia
using IterableTables, DataFrames
df = DataFrame(sol1)

# Result
17×5 DataFrames.DataFrame
│ Row │ timestamp │ value 1  │ value 2  │ value 3  │ value 4  │
├─────┼───────────┼──────────┼──────────┼──────────┼──────────┤
│ 1   │ 0.0       │ 0.110435 │ 0.569561 │ 0.918336 │ 0.508044 │
│ 2   │ 0.0625    │ 0.117406 │ 0.605515 │ 0.976306 │ 0.540114 │
│ 3   │ 0.125     │ 0.124817 │ 0.643738 │ 1.03794  │ 0.574208 │
│ 4   │ 0.1875    │ 0.132696 │ 0.684374 │ 1.10345  │ 0.610455 │
│ 5   │ 0.25      │ 0.141073 │ 0.727575 │ 1.17311  │ 0.64899  │
│ 6   │ 0.3125    │ 0.149978 │ 0.773503 │ 1.24716  │ 0.689958 │
│ 7   │ 0.375     │ 0.159445 │ 0.822331 │ 1.32589  │ 0.733511 │
│ 8   │ 0.4375    │ 0.16951  │ 0.87424  │ 1.40959  │ 0.779814 │
│ 9   │ 0.5       │ 0.18021  │ 0.929427 │ 1.49857  │ 0.82904  │
│ 10  │ 0.5625    │ 0.191586 │ 0.988097 │ 1.59316  │ 0.881373 │
│ 11  │ 0.625     │ 0.20368  │ 1.05047  │ 1.69373  │ 0.93701  │
│ 12  │ 0.6875    │ 0.216537 │ 1.11678  │ 1.80065  │ 0.996159 │
│ 13  │ 0.75      │ 0.230206 │ 1.18728  │ 1.91432  │ 1.05904  │
│ 14  │ 0.8125    │ 0.244738 │ 1.26222  │ 2.03516  │ 1.12589  │
│ 15  │ 0.875     │ 0.260187 │ 1.3419   │ 2.16363  │ 1.19697  │
│ 16  │ 0.9375    │ 0.276611 │ 1.42661  │ 2.30021  │ 1.27252  │
│ 17  │ 1.0       │ 0.294072 │ 1.51667  │ 2.44541  │ 1.35285  │
```

If a `ParameterizedFunction` is used, the output will use the variable names:

```julia
using ParameterizedFunctions

f = @ode_def begin
  dx = a*x - b*x*y
  dy = -3y + x*y
end a b

prob = ODEProblem(f,[1.0,1.0],(0.0,1.0),[1.5,1.0])
sol2 =solve(prob,Tsit5())

df = DataFrame(sol2)

7×3 DataFrames.DataFrame
│ Row │ timestamp │ x       │ y        │
├─────┼───────────┼─────────┼──────────┤
│ 1   │ 0.0       │ 1.0     │ 1.0      │
│ 2   │ 0.0776085 │ 1.04549 │ 0.857668 │
│ 3   │ 0.232645  │ 1.17587 │ 0.63946  │
│ 4   │ 0.429118  │ 1.41968 │ 0.456996 │
│ 5   │ 0.679082  │ 1.87672 │ 0.324733 │
│ 6   │ 0.944406  │ 2.58825 │ 0.263362 │
│ 7   │ 1.0       │ 2.77285 │ 0.25871  │
```

Additionally, this data can be saved to a CSV:

```julia
using CSV
CSV.write("out.csv",df)
```

For more information on using the IterableTables interface and other output
formats, see [IterableTables.jl](https://github.com/davidanthoff/IterableTables.jl).

## JLD2 and BSON.jl

JLD2.jl and BSON.jl will work with the full solution type if you bring the required functions
back into scope before loading. For eaxmple, if we save the solution:

```julia
using OrdinaryDiffEq, JLD2
f(u,p,t) = 1.01*u
u0=1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
@save "out.jld2" sol
```

then we can get the full solution type back, interpolations and all,
if we load the dependent functions first:

```julia
using JLD2
using OrdinaryDiffEq
f(u,p,t) = 1.01*u
JLD2.@load "out.jld2" sol
```

The example with BSON.jl is:

```julia
using OrdinaryDiffEq
f_2dlinear = (du,u,p,t) -> du.=1.01u
prob = ODEProblem(f_2dlinear,rand(2,2),(0.0,1.0))
sol1 =solve(prob,Euler();dt=1//2^(4))

using BSON
bson("test.bson",Dict(:sol1=>sol1))

# New session
using OrdinaryDiffEq
using BSON
BSON.load("test.bson")
```

If you load it without the DE function then for some algorithms the
interpolation may not work, and for all algorithms you'll need
at least a solver package or DiffEqBase.jl in scope in order for
the solution interface (plot recipes, array indexing, etc.) to
work. If none of these are put into scope, the solution type
will still load and hold all of the values (so `sol.u` and `sol.t`
will work), but none of the interface will be available.

## JLD

Don't use JLD. It's dead. Julia types can be saved via JLD.jl.
However, they cannot save types which have functions, which means that
the solution type is currently not compatible with JLD.

```julia
using JLD
JLD.save("out.jld","sol",sol)
```
