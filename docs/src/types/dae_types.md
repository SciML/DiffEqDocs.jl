# DAE Problems

```@docs
DAEProblem
DAEFunction
```

## Solution Type

```@docs
DAESolution
```

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/dae_premade_problems.jl).

To use a sample problem, such as `prob_dae_resrob`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.DAEProblemLibrary
# load problems
DAEProblemLibrary.importdaeproblems()
prob = DAEProblemLibrary.prob_dae_resrob
sol = solve(prob,IDA())
```

```@meta
CurrentModule = DAEProblemLibrary
```

```@docs
prob_dae_resrob
```
