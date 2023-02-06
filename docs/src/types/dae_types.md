# DAE Problems

```@docs
SciMLBase.DAEProblem
SciMLBase.DAEFunction
```

## Solution Type

```@docs
SciMLBase.DAESolution
```

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/SciML/DiffEqProblemLibrary.jl/blob/master/src/dae_premade_problems.jl).

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
