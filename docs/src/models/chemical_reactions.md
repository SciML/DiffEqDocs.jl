# Chemical Reactions

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) is a domain specific
 language (DSL) for the easy generation of models of chemical reaction systems,
 which can be used with SciML tooling to enable high performance simulation and
 analysis of chemical reaction networks. Catalyst generates `ReactionSystem`s,
 leveraging [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) to
 enable large-scale simulations through auto-vectorization and parallelism.
 `ReactionSystem`s can be used to generate ModelingToolkit-based models,
 allowing the easy simulation and parameter estimation of mass action ODE
 models, Chemical Langevin SDE models, stochastic chemical kinetics jump process
 models, and more. These generated models can be used with
 DifferentialEquations.jl solvers, but also with higher level SciML packages
 (e.g. for sensitivity analysis, parameter estimation, machine learning
 applications, etc). See the [Catalyst.jl
 documentation](https://catalyst.sciml.ai/) for more information.