# [Parameter Estimation and Bayesian Analysis](@id parameter_estimation)

Parameter estimation for differential equation models, also known as dynamic data analysis, 
is provided by the DiffEq suite. In this introduction, we briefly
present the relevant packages that facilitate parameter estimation, namely:

- [DiffEqFlux.jl](https://diffeqflux.sciml.ai/)
- [Turing.jl](https://turing.ml/)
- [DataDrivenDiffEq.jl](https://datadriven.sciml.ai/dev/)
- [DiffEqParamEstim.jl](https://diffeqparamestim.sciml.ai/dev/)
- [DiffEqBayes.jl](https://diffeqbayes.sciml.ai/dev/)

We also provide information regarding the respective strengths of these packages
so that you can easily decide which one suits your needs best.

### DiffEqFlux.jl

A very versatile and composable package, DiffEqFlux.jl allows for solving a
wide range of differential equations, for instance: stiff universal ODEs,
universal SDEs, universal PDEs, and other kinds of universal differential
equations. As regards probabilistic programming, DiffEqFlux.jl works in conjunction
with Turing.jl (see below). It is the most flexible and high performance
parameter estimation system.

### Turing.jl

In the context of differential equations and parameter estimation, Turing.jl
allows for a Bayesian estimation of differential equations (used in conjunction
with the high-level package DiffEqBayes.jl). For more examples on combining
Turing.jl with DiffEqBayes.jl, see the documentation below. It is important
to note that Turing.jl can also perform Bayesian estimation without relying on
DiffEqBayes.jl (for an example, consult [this](https://turing.ml/stable/tutorials/10-bayesian-differential-equations/) tutorial).

### DataDrivenDiffEq.jl

The distinguishing feature of this package is that its ultimate goal is to
identify the differential equation model that generated the input data.
Depending on the user's needs, the package can provide structural identification
of a given differential equation (output in a symbolic form) or structural
estimation (output as a function for prediction purposes).

### DiffEqParamEstim.jl

This package is for simplified parameter estimation. While not as flexible of a
system like DiffEqFlux.jl, it provides ready-made functions for doing standard
optmization procedures like L2 fitting and MAP estimates. Among other features,
it allows for the optimization of parameters in ODEs, stochastic problems, and
delay differential equations.

### DiffEqBayes.jl

As the name suggests, this package has been designed to provide the estimation
of differential equations parameters by means of Bayesian methods. It works in
conjunction with [Turing.jl](https://turing.ml/), 
[CmdStan.jl](https://github.com/StanJulia/CmdStan.jl), 
[DynamicHMC.jl](https://github.com/tpapp/DynamicHMC.jl), and 
[ApproxBayes.jl](https://github.com/marcjwilliams1/ApproxBayes.jl). While not
as flexible as direct usage of DiffEqFlux.jl or Turing.jl, DiffEqBayes.jl can
be an approachable interface for those not familiar with Bayesian estimation,
and provides a nice way to use Stan from pure Julia.
