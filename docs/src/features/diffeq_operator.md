# Matrix-Free Linear Operators with SciMLOperators.jl

SciML has the [SciMLOperators.jl](https://docs.sciml.ai/SciMLOperators/stable/) library for defining linear operators.
The ODE solvers will specialize on this property, for example using the linear operators automatically with Newton-Krylov
methods for matrix-free Newton-Krylov optimizations, but also allows for methods that require knowing that part of the
equation is linear, like exponential integrators. See the SciMLOperators.jl library for more information.
