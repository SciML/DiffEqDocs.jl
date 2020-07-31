# [Jump Problem and Jump Diffusion Solvers](@id jump_solve)

```julia
solve(prob::JumpProblem,alg;kwargs)
```

## Recommended Methods

A `JumpProblem(prob,aggregator,jumps...)` come in two forms. The first major
form is if it does not have a `RegularJump`. In this case, it can be solved with
any integrator on  `prob`. However, in the case of a pure `JumpProblem` (a
`JumpProblem` over a  `DiscreteProblem`), there are special algorithms
available.  The `SSAStepper()` is an efficient streamlined algorithm for running
the  `aggregator` version of the SSA for pure `ConstantRateJump` and/or
`MassActionJump` problems. However, it is not compatible with event handling. If
events are necessary, then `FunctionMap` does well.

If there is a `RegularJump`, then specific methods must be used. The current
recommended method is `TauLeaping` if you need adaptivity, events, etc. If you
just need the most barebones fixed time step leaping method, then `SimpleTauLeaping`
can have performance benefits.

## Special Methods for Pure Jump Problems

If you are using jumps with a differential equations, use the same methods
as in the case of the differential equation solving. However, the following
algorithms are optimized for pure jump problems.

### DiffEqJump.jl

- `SSAStepper`: a stepping algorithm for pure `ConstantRateJump` and/or
  `MassActionJump` `JumpProblem`s. Supports handling of `DiscreteCallback`
  and saving controls like `saveat`.

## RegularJump Compatible Methods

### StochasticDiffEq.jl

These methods support mixing with event handling, other jump types, and all of
the features of the normal differential equation solvers.

- `TauLeaping`: an adaptive tau-leaping algorithm with post-leap estimates.

### DiffEqJump.jl

- `SimpleTauLeaping`: a tau-leaping algorithm for pure `RegularJump` `JumpProblem`s.
  Requires a choice of `dt`.
- `RegularSSA`: a version of SSA for pure `RegularJump` `JumpProblem`s.

## Regular Jump Diffusion Compatible Methods

Regular jump diffusions are `JumpProblem`s where the internal problem is an `SDEProblem`
and the jump process has designed a regular jump.

### StochasticDiffEq.jl

- `EM`: Explicit Euler-Maruyama.
- `ImplicitEM`: Implicit Euler-Maruyama. See the SDE solvers page for more details.
