# BVP Problems

```@docs
SciMLBase.BVProblem
SciMLBase.SecondOrderBVProblem
```

## Solution Type

`BVProblem` solutions return an `ODESolution`. For more information, see the
[ODE problem definition page](@ref ode_prob) for the `ODESolution` docstring.

## Alias Specifier

```@docs
SciMLBase.BVPAliasSpecifier
```

## Example Problems

Example problems can be found in [DiffEqProblemLibrary.jl](https://github.com/SciML/DiffEqProblemLibrary.jl/blob/master/lib/BVProblemLibrary/src/BVProblemLibrary.jl).

To use a sample problem, such as `prob_bvp_linear_1`, you can do something like:

```julia
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.BVProblemLibrary, BoundaryValueDiffEq
# load problems
prob = BVProblemLibrary.prob_bvp_linear_1
sol = solve(prob, MIRK4(), dt = 0.05)
```

### Linear BVPs

```@meta
CurrentModule = BVProblemLibrary
```

```@docs
prob_bvp_linear_1
prob_bvp_linear_2
prob_bvp_linear_3
prob_bvp_linear_4
prob_bvp_linear_5
prob_bvp_linear_6
prob_bvp_linear_7
prob_bvp_linear_8
prob_bvp_linear_9
prob_bvp_linear_10
prob_bvp_linear_11
prob_bvp_linear_12
prob_bvp_linear_13
prob_bvp_linear_14
prob_bvp_linear_15
prob_bvp_linear_16
prob_bvp_linear_17
prob_bvp_linear_18
```

### Nonlinear BVPs

```@docs
prob_bvp_nonlinear_1
prob_bvp_nonlinear_2
prob_bvp_nonlinear_3
prob_bvp_nonlinear_4
prob_bvp_nonlinear_5
prob_bvp_nonlinear_6
prob_bvp_nonlinear_7
prob_bvp_nonlinear_8
prob_bvp_nonlinear_9
prob_bvp_nonlinear_10
prob_bvp_nonlinear_11
prob_bvp_nonlinear_12
prob_bvp_nonlinear_13
prob_bvp_nonlinear_14
prob_bvp_nonlinear_15
```

### Regular Nonlinear BVPs

```@docs
flat_moon
flat_earth
flat_earth_drag
measles
```
