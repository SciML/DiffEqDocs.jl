# [ODE Problems](@id ode_prob)

```@docs
SciMLBase.ODEProblem
SciMLBase.ODEFunction
```

## Solution Type

```@docs
SciMLBase.ODESolution
```

## Example Problems

Example problems can be found in [DiffEqProblemLibrary.jl](https://github.com/SciML/DiffEqProblemLibrary.jl/tree/master/src/ode).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.ODEProblemLibrary
# load problems
ODEProblemLibrary.importodeproblems()
prob = ODEProblemLibrary.prob_ode_linear
sol = solve(prob)
```

```@meta
CurrentModule = ODEProblemLibrary
```

```@docs
prob_ode_linear
prob_ode_2Dlinear
prob_ode_bigfloatlinear
prob_ode_bigfloat2Dlinear
prob_ode_large2Dlinear
prob_ode_2Dlinear_notinplace
prob_ode_lotkavolterra
prob_ode_fitzhughnagumo
prob_ode_threebody
prob_ode_pleiades
prob_ode_vanderpol
prob_ode_vanderpol_stiff
prob_ode_rober
prob_ode_rigidbody
prob_ode_hires
prob_ode_orego
prob_ode_pollution
prob_ode_nonlinchem
prob_ode_brusselator_1d
prob_ode_brusselator_2d
prob_ode_filament
prob_ode_thomas
prob_ode_lorenz
prob_ode_aizawa
prob_ode_dadras
prob_ode_chen
prob_ode_rossler
prob_ode_rabinovich_fabrikant
prob_ode_sprott
prob_ode_hindmarsh_rose
```
