# [Local Sensitivity Analysis (Automatic Differentiation)](@id sensitivity)

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

There are three types of sensitivity analysis. Local forward sensitivity
analysis directly gives the gradient of the solution with respect to each
parameter along the time series. The computational cost scales like `N*M`,
where `N` is the number of states and `M` is the number of parameters. While
this gives all of the information, it can be expensive for models with large
numbers of parameters. Local adjoint sensitivity analysis solves directly for
the gradient of some functional of the solution, such as a cost function or
energy functional, in a manner that is cheaper when the number of parameters is
large. Global Sensitivity Analysis methods are meant to be used for exploring the
sensitivity over a larger domain without calculating derivatives and are covered
on a different page.

## Installation and Usage

This functionality does not come standard with DifferentialEquations.jl.
To use this functionality, you must install SciMLSensitivity.jl:

```julia
]add SciMLSensitivity
using SciMLSensitivity
```

For complete information on using the sensitivity analyis features, please
[consult the SciMLSensitivity.jl documentation](https://sensitivity.sciml.ai/dev)
