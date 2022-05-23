# [Dynamical, Hamiltonian and 2nd Order ODE Problems](@id dynamical_prob)

```@docs
DynamicalODEProblem
SecondOrderODEProblem
DynamicalODEFunction
```

## Solution Type

Dynamical ODE solutions return an `ODESolution`. For more information, see the
[ODE problem definition page](@ref ode_prob) for the `ODESolution` docstring.

## Hamiltonian Problems

`HamiltonianProblem`s are provided by DiffEqPhysics.jl and provide an easy way
to define equations of motion from the corresponding Hamiltonian. To define a
`HamiltonianProblem` one only needs to specify the Hamiltonian:

```math
H(p,q)
```

and autodifferentiation (via ForwardDiff.jl) will create the appropriate
equations.

### Constructors

```julia
HamiltonianProblem{T}(H,p0,q0,tspan,param=nothing;kwargs...)
```

### Fields

* `H`: The Hamiltonian `H(p,q,params)` which returns a scalar.
* `p0`: The initial momentums.
* `q0`: The initial positions.
* `tspan`: The timespan for the problem.
* `param`: Defaults to `nothing`. `param` will be passed to `H`'s `params`. 
