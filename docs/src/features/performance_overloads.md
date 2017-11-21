# Performance Overloads

The DiffEq ecosystem provides an extensive interface for declaring extra functions
associated with the differential equation's data. In traditional libraries there
is usually only one option: the Jacobian. However, we allow for a large array
of pre-computed functions to speed up the calculations. This is offered via function
overloading (or overloaded types) and allows for these extra features to be
offered without cluttering the problem interface.

## Declaring Explicit Jacobians

The most standard case, declaring a function for a Jacobian is done by overloading
the function `f(t,u,du)` with an in-place updating function for the Jacobian:
`f(Val{:jac},t,u,J)` where the value type is used for dispatch. For example,
take the LotkaVolterra model:

```julia
function f(t,u,du)
  du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
end
```

To declare the Jacobian we simply add the dispatch:

```julia
function f(::Type{Val{:jac}},t,u,J)
  J[1,1] = p.a - p.b * u[2]
  J[1,2] = -(p.b) * u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  nothing
end
```

Note that this can also be done by generating a call-overloaded type. Indeed, this
is what ParameterizedFunctions.jl does, so see [its README](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl).

## Declaring Explicit Jacobians for DAEs

For fully implicit ODEs (`DAEProblem`s), a slightly different Jacobian function
is necessary. For the DAE

```math
G(t,u,du) = res
```

The Jacobian should be given in the form `dG/du + gamma*dG/d(du)` where `gamma`
is given by the solver. This means that the signature is:

```julia
f(::Type{Val{:jac}},t,u,du,gamma,J)
```

For example, for the equation

```julia
function testjac(t,u,du,res)
  res[1] = du[1] - 2.0 * u[1] + 1.2 * u[1]*u[2]
  res[2] = du[2] -3 * u[2] - u[1]*u[2]
end
```

we would define the Jacobian as:

```julia
function testjac(::Type{Val{:jac}},t,u,du,gamma,J)
  J[1,1] = gamma - 2.0 + 1.2 * u[2]
  J[1,2] = 1.2 * u[1]
  J[2,1] = - 1 * u[2]
  J[2,2] = gamma - 3 - u[1]
  nothing
end
```

## Other Available Functions

The full interface available to the solvers is as follows:

```julia
f(t,u,du) # Call the function
f(Val{:analytic},t,u,du) # The analytical solution. Used in testing
f(t,u,params,du) # Call the function to calculate with parameters params (vector)
f(Val{:tgrad},t,u,J) # Call the explicit t-gradient function
f(Val{:paramjac},t,u,params,J) # Call the explicit parameter Jacobian function
f(Val{:jac},t,u,J) # Call the explicit Jacobian function
f(Val{:invjac},t,u,iJ) # Call the explicit Inverse Jacobian function
f(Val{:invW},t,u,γ,iW) # Call the explicit inverse Rosenbrock-W function (M - γJ)^(-1)
f(Val{:invW_t},t,u,γ,iW) # Call the explicit transformed inverse Rosenbrock-W function (M/γ - J)^(-1)
```

Overloads which require parameters should subtype `ParameterizedFunction`. These
are all in-place functions which write into the last variable. See solver documentation
specifics to know which optimizations the algorithms can use.

## Symbolically Calculating the Functions

ParameterizedFunctions.jl automatically calculates as many of these functions as
possible and generates the overloads using SymEngine. Thus, for best performance
with the least work, it is suggested one use ParameterizedFunctions.jl.
