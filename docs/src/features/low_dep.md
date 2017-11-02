# Low Dependency Usage

DifferentialEquations.jl is a large library containing the functionality of
many different solver and addon packages. However in many cases you may want
to cut down on the size of the dependency and only use the parts of the
the library which are essential to your application. This is possible
due to JuliaDiffEq's modular package structure.

## Common Example: Using only OrdinaryDiffEq.jl

One common example is using only the ODE solvers OrdinaryDiffEq.jl. The solvers all
reexport DiffEqBase.jl (which holds the problem and solution types) and so
OrdinaryDiffEq.jl is all that's needed. Thus replacing

```julia
using DifferentialEquations
```

with

```julia
using OrdinaryDiffEq
```

will work if these are the only features you are using.

## Generalizing the Idea

In general, you will always need DiffEqBase.jl, since it defines all of the
fundamental types, but the solvers will automatically reexport it. 
For solvers, you typically only need that solver package.
So DiffEqBase+Sundials, DiffEqBase+LSODA, etc. will get you the common interface
with that specific solver setup. DiffEqBase.jl is a very lightweight dependency,
so there is no issue here! For PDEs, you normally need DiffEqBase+DiffEqPDEBase
in addition to the solver package.

For the addon packages, you will normally need DiffEqBase, the solver package
you choose, and the addon package. So for example, for parameter estimation you
would likely want DiffEqBase+OrdinaryDiffEq+DiffEqParamEstim. If you arne't sure
which package a specific command is from, they using `@which`. For example, from
the parameter estimation docs we have:

```julia
using DifferentialEquations
f = @ode_def_nohes LotkaVolterraTest begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=1.0 c=3.0 d=1.0

u0 = [1.0;1.0]
tspan = (0.0,10.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5())
t = collect(linspace(0,10,200))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
using RecursiveArrayTools
data = convert(Array,randomized)
cost_function = build_loss_objective(prob,t,data,Tsit5(),maxiters=10000)
```

If we wanted to know where `build_loss_objective` came from, we can do:

```julia
@which build_loss_objective(prob,t,data,Tsit5(),maxiters=10000)

(::DiffEqParamEstim.#kw##build_loss_objective)(::Array{Any,1}, ::DiffEqParamEstim.#build_loss_objective, prob::DiffEqBase.DEProblem, t, data, alg)
```

This says it's in the DiffEqParamEstim.jl package. Thus in this case, we could have
done

```julia
using OrdinaryDiffEq, DiffEqParamEstim
```

instead of the full `using DifferentialEquations`. Note that due to the way
Julia dependencies work, any internal function in the package will work. The only
dependencies you need to explicitly `using` are the functions you are specifically
calling. Thus this method can be used to determine all of the DiffEq packages
you are using.
