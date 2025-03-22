# Matrix-Free Linear Operators and Specializations on Linearity

SciML has the [SciMLOpereators.jl](https://docs.sciml.ai/SciMLOperators/stable/) library for defining linear operators.
The ODE solvers will specialize on this property, for example using the lienar operators automatically with Newton-Krylov
methods for matrix-free Newton-Krylov optimzations, but also allows for methods that requires knowing that part of the
equation is linear, like exponential integrators. See the SciMLOperators.jl library for more information.
