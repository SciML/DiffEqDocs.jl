# SDE Problems

```@docs
SciMLBase.SDEProblem
SciMLBase.SDEFunction
```

## Solution Type

`SDEProblem` solutions return an `RODESolution`. For more information, see the
[RODE problem definition page](@ref rode_problem) for the `RODESolution` docstring.

## Alias Specifier

```@docs
SciMLBase.SDEAliasSpecifier
```

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/SciML/DiffEqProblemLibrary.jl/blob/master/lib/SDEProblemLibrary/src/SDEProblemLibrary.jl).

To use a sample problem, such as `prob_sde_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
import DiffEqProblemLibrary.SDEProblemLibrary
import DifferentialEquations as DE
# load problems
SDEProblemLibrary.importsdeproblems()
prob = SDEProblemLibrary.prob_sde_linear
sol = DE.solve(prob)
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
