# [Jump Problem Solvers](@id jump_solve)

`solve(prob::JumpProblem,alg;kwargs)`

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
recommended method is `SimpleTauLeaping`.

## Special Methods for Pure Jump Problems

If you are using jumps with a differential equations, use the same methods
as in the case of the differential equation solving. However, the following
algorithms are optimized for pure jump problems.

### DiffEqJump.jl

- `SSAStepper`: a stepping algorithm for pure `ConstantRateJump` and/or
  `MassActionJump` `JumpProblem`s. Does not support event handling, but does
  support saving controls like `saveat`.

## RegularJump Compatible Methods

### DiffEqJump.jl

- `SimpleTauLeaping`: a tau-leaping algorithm for pure `RegularJump` `JumpProblem`s.
  Requires a choice of `dt`.
- `RegularSSA`: a version of SSA for pure `RegularJump` `JumpProblem`s.
