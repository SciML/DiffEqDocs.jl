# [Reduced Compile Time, Optimizing Runtime, and Low Dependency Usage](@id low_dep)

While DifferentialEquations.jl's defaults strike a balance between runtime
and compile time performance, users should be aware that there are many
options for controlling the dependency sizing and the compilation caching
behavior to further refine this trade-off. The following methods are
available for such controls.

## Controlling Function Specialization and Precompilation

By default, DifferentialEquations.jl solvers make use of function wrapping
techniques in order to fully precompile the solvers and thus decrease the
compile time. However, in some cases you may wish to control this behavior,
pushing more towards faster runtimes or faster compile times. This can be
done by using the `specialization` arguments of the `AbstractDEProblem`
constructors.

For example, with the `ODEProblem` we have `ODEProblem{iip,specialize}(...)`.
This second type parameter controls the specialization level with the
following choices:

  - `SciMLBase.AutoSpecialize`: the default. Uses a late wrapping scheme to
    hit a balance between runtime and compile time.
  - `SciMLBase.NoSpecialize`: this will never specialize on the constituent
    functions, having the least compile time but the highest runtime.
  - `SciMLBase.FullSpecialize`: this will fully re-specialize the solver
    on the given ODE, achieving the fastest runtimes while increasing the
    compile times. This is what is recommended when benchmarking and when
    running long computations, such as in optimization loops.

For more information on the specialization levels, please see
[the SciMLBase documentation on specialization levels](https://docs.sciml.ai/SciMLBase/stable/interfaces/Problems/#Specialization-Levels).

DifferentialEquations.jl and its ODE package OrdinaryDiffEq.jl precompile
some standard problem types and solvers. The problem types include the
three specialization levels described above and the default setting.
The solvers include some

  - standard solvers for non-stiff problems such as `Tsit5()`
  - standard solvers for stiff problems such as `Rosenbrock23()`
  - standard solvers with stiffness detection such as `AutoTsit5(Rosenbrock23())`
  - low-storage methods for conservation laws such as `SSPRK43()`
    (precompilation disabled by default)

To adapt the amount of precompilation, you can use
[Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl).
For example, to turn off precompilation for non-default problem types
(specialization levels) and all stiff/implicit/low-storage solvers,
you can execute the following code in your active project.

```
using Preferences, UUIDs
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileNonStiff" => true)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileStiff" => false)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileAutoSwitch" => false)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileLowStorage" => false)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileDefaultSpecialize" => true)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileAutoSpecialize" => false)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileFunctionWrapperSpecialize" => false)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileNoSpecialize" => false)
```

This will create a `LocalPreferences.toml` file next to the currently active
`Project.toml` file.

## Decreasing Dependency Size by Direct Dependence on Specific Solvers

DifferentialEquations.jl is a large library containing the functionality of
many different solver and add-on packages. However, often you may want
to cut down on the size of the dependency and only use the parts of
the library which are essential to your application. This is possible
due to SciML's modular package structure.

### Common Example: Using only OrdinaryDiffEq.jl

One common example is using only the ODE solvers OrdinaryDiffEq.jl. The solvers all
reexport SciMLBase.jl (which holds the problem and solution types) and so
OrdinaryDiffEq.jl is all that's needed. Thus replacing

```julia
using DifferentialEquations
```

with

```julia
#Add the OrdinaryDiffEq Package first!
#using Pkg; Pkg.add("OrdinaryDiffEq")
using OrdinaryDiffEq
```

will work if these are the only features you are using.

### Generalizing the Idea

In general, you will always need SciMLBase.jl, since it defines all of the
fundamental types, but the solvers will automatically reexport it.
For solvers, you typically only require that solver package.
So SciMLBase+Sundials, SciMLBase+LSODA, etc. will get you the common interface
with that specific solver setup. SciMLBase.jl is a very lightweight dependency,
so there is no issue here!

For the add-on packages, you will normally need SciMLBase, the solver package
you choose, and the add-on package. So for example, for predefined callbacks you
would likely want SciMLBase+OrdinaryDiffEq+DiffEqCallbacks. If you aren't sure
which package a specific command is from, then use `@which`. For example, from
the callback docs we have:

```@example low_dep_1
using DifferentialEquations
function fitz(du, u, p, t)
    V, R = u
    a, b, c = p
    du[1] = c * (V - V^3 / 3 + R)
    du[2] = -(1 / c) * (V - a - b * R)
end
u0 = [-1.0; 1.0]
tspan = (0.0, 20.0)
p = (0.2, 0.2, 3.0)
prob = ODEProblem(fitz, u0, tspan, p)
cb = ProbIntsUncertainty(0.2, 1)
ensemble_prob = EnsembleProblem(prob)
sim = solve(ensemble_prob, Euler(), trajectories = 100, callback = cb, dt = 1 / 10)
```

If we wanted to know where `ProbIntsUncertainty(0.2,1)` came from, we can do:

```@example low_dep_1
using InteractiveUtils # hide
@which ProbIntsUncertainty(0.2, 1)
```

This says it's in the DiffEqCallbacks.jl package. Thus in this case, we could have
done

```@example low_dep_2
using OrdinaryDiffEq, DiffEqCallbacks
function fitz(du, u, p, t)
    V, R = u
    a, b, c = p
    du[1] = c * (V - V^3 / 3 + R)
    du[2] = -(1 / c) * (V - a - b * R)
end
u0 = [-1.0; 1.0]
tspan = (0.0, 20.0)
p = (0.2, 0.2, 3.0)
prob = ODEProblem(fitz, u0, tspan, p)
cb = ProbIntsUncertainty(0.2, 1)
ensemble_prob = EnsembleProblem(prob)
sim = solve(ensemble_prob, Euler(), trajectories = 100, callback = cb, dt = 1 / 10)
```

instead of the full `using DifferentialEquations`. Note that due to the way
Julia dependencies work, any internal function in the package will work. The only
dependencies you need to explicitly `using` are the functions you are specifically
calling. Thus, this method can be used to determine all of the DiffEq packages
you are using.
