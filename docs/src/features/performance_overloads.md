# DiffEqFunctions (Jacobians, Gradients, etc.) and Jacobian Types

The DiffEq ecosystem provides an extensive interface for declaring extra functions
associated with the differential equation's data. In traditional libraries there
is usually only one option: the Jacobian. However, we allow for a large array
of pre-computed functions to speed up the calculations. This is offered via the
`DiffEqFunction` types which can be passed to the problems.

## Function Type Definitions

### Function Choice Definitions

The full interface available to the solvers is as follows:

- `jac`: The Jacobian of the differential equation with respect to the state
  variable `u` at a time `t` with parameters `p`.
- `tgrad`: The gradient of the differential equation with respect to `t` at state
  `u` with parameters `p`.
- `paramjac`: The Jacobian of the differential equation with respect to `p` at
  state `u` at time `t`.
- `analytic`: Defines an analytical solution using `u0` at time `t` with `p`
  which will cause the solvers to return errors. Used for testing.
- `Wfact`: The LU-factorization of `M - gamma*J` where `J` is the `jac`.
- `Wfact_t`: The LU-factorization of `M/gamma - J` where `J` is the `jac`.
- `ggprime`: See the definition in the SDEProblem page.
- `syms`: Allows you to name your variables for automatic names in plots and
  other output.

### ODEFunction

```julia
function ODEFunction{iip,recompile}(f;
                 mass_matrix=I,
                 analytic=nothing, # (u0,p,t)
                 tgrad=nothing, # (dT,u,p,t) or (u,p,t)
                 jac=nothing, # (J,u,p,t) or (u,p,t)
                 jac_prototype=nothing, # Type for the Jacobian
                 Wfact=nothing, # (iW,u,p,gamma,t) or (u,p,gamma,t)
                 Wfact_t=nothing, # (iW,u,p,gamma,t) or (u,p,gamma,t)
                 paramjac = nothing, # (pJ,u,p,t) or (u,p,t)
                 colorvec = nothing,
                 syms = nothing) # collection of names for variables
```

### DynamicalODEFunction

```julia
DynamicalODEFunction{iip,recompile}(f1, # (du,u,v,p,t) or (u,v,p,t)
                                    f2; # (du,u,v,p,t) or (u,v,p,t)
                                    mass_matrix=(I,I), # Mass matrix for each partition
                                    analytic=nothing)
```

### SplitFunction

```julia
SplitFunction{iip,recompile}(f1, # ODEFunction
                        f2; # ODEFunction
                        mass_matrix=I,
                        _func_cache=nothing, # This is a cache used in f = f1+f2
                        analytic=nothing)
```

### SDEFunction

```julia
function SDEFunction{iip,recompile}(f,g;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 colorvec = nothing,
                 syms = nothing)
```

### SplitSDEFunction

```julia
SplitSDEFunction{iip,recompile}(f1, # ODEFunction
                           f2, # ODEFunction
                           g;
                           mass_matrix=I,
                           _func_cache=nothing,
                           analytic=nothing)
```

### RODEFunction

```julia
function RODEFunction{iip,recompile}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 colorvec = nothing,
                 syms = nothing)
```

### DAEFunction

```julia
function DAEFunction{iip,recompile}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing, # (J,du,u,p,gamma,t) or (du,u,p,gamma,t)
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing)
```

Note that the Jacobian of a DAE is defined as `gamma*dG/d(du) + dG/du ` where
`gamma` is given by the solver.

### DDEFunction

```julia
function DDEFunction{iip,recompile}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 colorvec = nothing,
                 syms = nothing)
```

## Inplace Specification and No-Recompile Mode

Each DiffEqFunction type can be called with an "is inplace" (iip) choice.

```julia
ODEFunction(f)
ODEFunction{iip}(f)
```

which is a boolean for whether the function is in the inplace form (mutating to
change the first value). This is automatically determined using the methods table
but note that for full type-inferrability of the `DEProblem` this iip-ness should
be specified.

Additionally, the functions are fully specialized to reduce the runtimes. If one
would instead like to not specialize on the functions to reduce compile time,
then one can set `recompile` to false.

```julia
ODEFunction{iip,false}(f)
```

This makes the ODE solver compilation independent of the function and so changing
the function will not cause recompilation. One can change the default value
by changing the `const RECOMPILE_BY_DEFAULT = true` to false in the DiffEqBase.jl
source code.

## Specifying Jacobian Types

The `jac` field of an inplace style `DiffEqFunction` has the signature `jac(J,u,p,t)`,
which updates the jacobian `J` in-place. The intended type for `J` can sometimes be
inferred (e.g. when it is just a dense `Matrix`), but not in general. To supply the
type information, you can provide a `jac_prototype` in the function's constructor.

The following example creates an inplace `ODEFunction` whose jacobian is a `Diagonal`:

```julia
using LinearAlgebra
f = (du,u,p,t) -> du .= t .* u
jac = (J,u,p,t) -> (J[1,1] = t; J[2,2] = t; J)
jp = Diagonal(zeros(2))
fun = ODEFunction(f; jac=jac, jac_prototype=jp)
```

Note that the integrators will always make a deep copy of `fun.jac_prototype`, so
there's no worry of aliasing.

In general the jacobian prototype can be anything that has `mul!` defined, in
particular sparse matrices or custom lazy types that support `mul!`. A special case
is when the `jac_prototype` is a `AbstractDiffEqLinearOperator`, in which case you
do not need to supply `jac` as it is automatically set to `update_coefficients!`.
Refer to the [DiffEqOperators](../diffeq_operator) section for more information
on setting up time/parameter dependent operators.

## Examples

### Declaring Explicit Jacobians for ODEs

The most standard case, declaring a function for a Jacobian is done by overloading
the function `f(du,u,p,t)` with an in-place updating function for the Jacobian:
`f_jac(J,u,p,t)` where the value type is used for dispatch. For example,
take the LotkaVolterra model:

```julia
function f(du,u,p,t)
  du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
end
```

To declare the Jacobian we simply add the dispatch:

```julia
function f_jac(J,u,p,t)
  J[1,1] = 2.0 - 1.2 * u[2]
  J[1,2] = -1.2 * u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  nothing
end
```

Then we can supply the Jacobian with our ODE as:

```julia
ff = ODEFunction(f;jac=f_jac)
```

and use this in an `ODEProblem`:

```julia
prob = ODEProblem(ff,ones(2),(0.0,10.0))
```

### Declaring Explicit Jacobians for DAEs

For fully implicit ODEs (`DAEProblem`s), a slightly different Jacobian function
is necessary. For the DAE

```math
G(du,u,p,t) = res
```

The Jacobian should be given in the form `gamma*dG/d(du) + dG/du ` where `gamma`
is given by the solver. This means that the signature is:

```julia
f(J,du,u,p,gamma,t)
```

For example, for the equation

```julia
function testjac(res,du,u,p,t)
  res[1] = du[1] - 2.0 * u[1] + 1.2 * u[1]*u[2]
  res[2] = du[2] -3 * u[2] - u[1]*u[2]
end
```

we would define the Jacobian as:

```julia
function testjac(J,du,u,p,gamma,t)
  J[1,1] = gamma - 2.0 + 1.2 * u[2]
  J[1,2] = 1.2 * u[1]
  J[2,1] = - 1 * u[2]
  J[2,2] = gamma - 3 - u[1]
  nothing
end
```

## Symbolically Calculating the Functions

See the `modelingtoolkitize` function from
[ModelingToolkit.jl](https://github.com/JuliaDiffEq/ModelingToolkit.jl) for
automatically symbolically calculating the Jacobian for numerically-defined
functions.
