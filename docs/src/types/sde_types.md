# SDE Problems

```@docs
SDEProblem
```

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/sde_premade_problems.jl).

To use a sample problem, such as `prob_sde_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.SDEProblemLibrary
# load problems
SDEProblemLibrary.importsdeproblems()
prob = SDEProblemLibrary.prob_sde_linear
sol = solve(prob)
```

```@meta
CurrentModule = SDEProblemLibrary
```

```@docs
prob_sde_linear
prob_sde_2Dlinear
prob_sde_wave
prob_sde_lorenz
prob_sde_cubic
prob_sde_additive
prob_sde_additivesystem
prob_sde_nltest
oval2ModelExample
prob_sde_stiffquadstrat
prob_sde_stiffquadito
generate_stiff_stoch_heat
prob_sde_bistable
prob_sde_bruss
prob_sde_oscilreact
```
