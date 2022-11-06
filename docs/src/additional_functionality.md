# Additional Functionality
This page gives a short overview of other tools that have been integrated with DifferentialEquations.jl.

## Domain Modeling Tools
Many domain specific modeling tools have been build atop of DifferentialEquations.jl,
you can find them [here](https://docs.sciml.ai/Overview/stable/highlevels/modeling_languages/).

## Parameter Estimation and Bayesian Analysis

Parameter estimation for differential equation models, also known as dynamic data analysis
and as inverse problems, is tackled by
[various packages in the Sci-ML ecosystem.](https://docs.sciml.ai/Overview/stable/highlevels/inverse_problems/).

## Neural Networks

Use DifferentialEquations.jl with various
[neural network packages](https://docs.sciml.ai/Overview/stable/highlevels/function_approximation/#Third-Party-Libraries-to-Note),
such as Lux.jl.

## Bifurcation Analysis

Bifurcation analysis on DifferentialEquations.jl types can be performed by,
[BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl),
BifurcationKit has integration with `ODEProblem` for some functionality (like computing periodic orbits via shooting).
If `oprob` is an `ODEProblem`, one can also just pass `oprob.f.f` and `oprob.f.jac` to BifurcationKit methods as needed.
Bifurcations.jl can directly generate a `BifurcationProblem` from an
`ODEProblem`. 

## Uncertainty Quantification

There's always uncertainty in our models. For various tools to quantify
this uncertainty, see 
[the DiffEqCallbacks documentation](https://docs.sciml.ai/Overview/stable/highlevels/uncertainty_quantification/).

## Global Sensitivity Analysis

Global Sensitivity Analysis (GSA) methods are used to quantify the uncertainty in
output of a model w.r.t. the parameters, their individual contributions, or the
contribution of their interactions. The GSA interface allows for utilizing batched
functions for parallel computation of GSA quantities.

For more information, see the documentation on [GlobalSensitivity.jl](https://docs.sciml.ai/GlobalSensitivity/stable/).

## Local Sensitivity Analysis (Automatic Differentiation)

Sensitivity analysis, or automatic differentiation of the solver, is provided
by the DiffEq suite. The model sensitivities are the derivatives of the
solution ``u(t)`` with respect to the parameters. Specifically, the local
sensitivity of the solution to a parameter is defined by how much the solution
would change by changes in the parameter, i.e. the sensitivity of the ith
independent variable to the jth parameter is
`` \frac{\partial u_i}{\partial p_{j}}``.

Sensitivity analysis serves two major purposes. On one hand, the sensitivities
are diagnostics of the model which are useful for understand how
it will change in accordance to changes in the parameters. But another use is
simply because in many cases these derivatives are useful. Sensitivity analysis
provides a cheap way to calculate the gradient of the solution which can be
used in parameter estimation and other optimization tasks.

There are two types of sensitivity analysis. Local forward sensitivity
analysis directly gives the gradient of the solution with respect to each
parameter along the time series. The computational cost scales like `N*M`,
where `N` is the number of states and `M` is the number of parameters. While
this gives all of the information, it can be expensive for models with large
numbers of parameters. Local adjoint sensitivity analysis solves directly for
the gradient of some functional of the solution, such as a cost function or
energy functional, in a manner that is cheaper when the number of parameters is
large.

This functionality does not come standard with DifferentialEquations.jl.
To use this functionality, you must install SciMLSensitivity.jl:

```julia
using Pkg
Pkg.add("SciMLSensitivity")
using SciMLSensitivity
```

For complete information on using the sensitivity analyis features, please
[consult the SciMLSensitivity.jl documentation](https://docs.sciml.ai/SciMLSensitivity/stable/)
