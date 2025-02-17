# DAE Problems

```@docs
SciMLBase.DAEProblem
SciMLBase.DAEFunction
```

## Solution Type

```@docs
SciMLBase.DAESolution
```

## Alias Specifier

```@docs
SciMLBase.DAEAliasSpecifier
```

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/SciML/DiffEqProblemLibrary.jl/tree/master/lib/DAEProblemLibrary).

To use a sample problem, such as `prob_dae_resrob`, you can do something like:

```julia
using DiffEqProblemLibrary.DAEProblemLibrary, Sundials
prob = DAEProblemLibrary.prob_dae_resrob
sol = solve(prob, IDA())
```

```@meta
CurrentModule = DAEProblemLibrary
```

```@docs
prob_dae_resrob
```
