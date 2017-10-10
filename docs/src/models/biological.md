# Chemical Reaction Models

The biological models functionality is provided by DiffEqBiological.jl and helps
the user build discrete stochastic and differential equation based systems
biological models. These tools allow one to define the models at a high level
by specifying reactions and rate constants, and the creation of the actual problems
is then handled by the modeling package.

## The Reaction Type

The basic type for BiologicalModels.jl is the reaction type. Its constructor is:

```julia
Reaction(rate_constant,reactants,stoichiometry)
```

`rate_constant` is the rate constant for the reaction. `reactants` is a list of
reactants for the reaction. For example, `reactants=(2,3)` means that the reaction
rate is `rate_constant*u[2]*u[3]`. `stoichiometry` is then the stoichiometry for
the result. It is a list of tuples of changes to apply when the reaction takes place.
Each tuple `(i,j)` means "modify reactiant i by amount j". For example, the tuple
`(2,-1)` means "decrease reactant 2 by 1".

## Note About Rate Dependency

Note that currently, the reactions are used to build `ConstantRateJump`s. This means
that the solver requires that the rates are constant between jumps in order to
achieve full accuracy. The rates for the `ConstantRateJump` may depend on each other,
but they may not depend on the differential equation themselves.

## Variable Rate Reactions

`VariableRateReaction` are allowed to have their rates change continuously, depending
on time or values related to a differential equation. The constructor is:

```julia
function VariableRateReaction(rate_constant,reactants,stoichiometry;
                              idxs = nothing,
                              rootfind=true,
                              interp_points=10,
                              abstol=1e-12,reltol=0)
```

The additional keyword arguments are for controlling the associated `ContinuousCallback`
used to handle `VariableRateReaction`s in simulations.

### Example Reaction

An example reaction is:

```julia
r1 = Reaction(1e-4,(1,2),((1,-1),(2,1)))
```

Here, the `rate_constant` is `1e-4`. The reactants are components 1 and 2, meaning
the reaction rate is calculated by `rate_constant*u[1]*u[2]`. The `stoichiometry`
does two things. First, the `(1,-1)` means that, when the reaction occurs, we
decrease `u[1]` by 1. Secondly, the `(2,1)` means we increase `u[2]` by 1. Thus
this reaction is a reaction where chemical 1 changes into chemical 2, and it is
enhanced by chemical 2 itself.

## GillespieProblem

These reactions can be added to a differential equation (or discrete) problem
using the `GillespieProblem`. This is simply a constructor which interprets the
reactions as jumps, and builds the associated `JumpProblem`. Thus its constructor
is the same:

```julia
GillespieProblem(prob,aggregator::AbstractAggregatorAlgorithm,rs::AbstractReaction...;kwargs...)
```

This is the exact same constructor as the `JumpProblem`, except now we pass reactions
(or `VariableRateReaction`s, or a `ReactionSet`) instead of jumps. Thus for more
information, see the description of the `JumpProblem`.

## The Reaction DSL

The `@reaction_network` DSL allows you to define reaction networks in a more
scientific format. Each line is given as `parameter reactants --> products`.

### Example: Birth-Death Process

```julia
rs = @reaction_network begin
  2.0, X --> 2X
  1.0, X --> 0
  0.5, 0 --> X
end
prob = DiscreteProblem([5], (0.0, 4.0))
jump_prob = GillespieProblem(prob, Direct(), rs)
sol = solve(jump_prob, Discrete())
```

### Example: Michaelis-Menten Enzyme Kinetics

```julia
rs = @reaction_network begin
  0.00166, S + E --> SE
  0.0001,  SE --> S + E
  0.1,     SE --> P + E
end
# S = 301, E = 100, SE = 0, P = 0
prob = DiscreteProblem([301, 100, 0, 0], (0.0, 100.0))
jump_prob = GillespieProblem(prob, Direct(), rs)
sol = solve(jump_prob, Discrete())
```
