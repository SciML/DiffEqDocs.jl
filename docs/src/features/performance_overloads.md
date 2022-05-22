# [Jacobians, Gradients, etc.](@id performance_overloads)

The DiffEq ecosystem provides an extensive interface for declaring extra functions
associated with the differential equation's data. In traditional libraries there
is usually only one option: the Jacobian. However, we allow for a large array
of pre-computed functions to speed up the calculations. This is offered via the
`DiffEqFunction` types which can be passed to the problems.

## Built-In Jacobian Options

!!! warning

    This subsection on Jacobian options only applies to the Julia-based solvers
    (OrdinaryDiffEq.jl, StochasticDiffEq.jl, etc.). Wrappers of traditional
    C/Fortran codes, such as Sundials.jl, always default to finite differencing
    using the in-solver code.

All applicable stiff differential equation solvers in the Julia ecosystem
(OrdinaryDiffEq.jl, StochasticDiffEq.jl, DelayDiffEq.jl, etc.) take the following
arguments for handling the automatic Jacobian construction with the following defaults:

- `chunk_size`: The chunk size used with ForwardDiff.jl. Defaults to `Val{0}()`
  and thus uses the internal ForwardDiff.jl algorithm for the choice.
- `autodiff`: Specifies whether to use automatic differentiation via
  [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) or finite
  differencing via [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl).
  Defaults to `Val{true}()` for automatic differentation.
- `standardtag`: Specifies whether to use package-specific tags instead of the
  ForwardDiff default function-specific tags. For more information see
  [this blog post](https://www.stochasticlifestyle.com/improved-forwarddiff-jl-stacktraces-with-package-tags/).
  Defaults to `Val{true}()`.
- `concrete_jac`: Specifies whether a Jacobian should be constructed. Defalts to
  `nothing`, which means it will be chosen true/false depending on circumstances
  of the solver, such as whether a Krylov subspace method is used for `linsolve`.
- `diff_type`: The type of differentiation used in FiniteDiff.jl if `autodiff=false`.
  Defalts to `Val{:forward}`, with altnerative choices of `Val{:central}` and
  `Val{:complex}`. 

## Function Type Definitions

```@docs
ODEFunction
SDEFunction
DAEFunction
DiscreteFunction
DDEFunction
SDDEFunction
DynamicalODEFunction
DynamicalSDEFunction
DynamicalDDEFunction
SplitFunction
SplitSDEFunction
```