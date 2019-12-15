# [ParameterizedFunctions](@id paremeterized_functions)

## Installation

This functionality does not come standard with DifferentialEquations.jl.
To use this functionality, you must install ParameterizedFunctions.jl:

```julia
]add ParameterizedFunctions
using ParameterizedFunctions
```

## Function Definition Macros

DifferentialEquations.jl provides a set of macros for more easily and legibly
defining your differential equations. It exploits the standard notation for
mathematically writing differential equations and the notation for "punching
differential equations into the computer"; effectively doing the translation
step for you. This is best shown by an example. Say we want to solve the
[ROBER model](http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Single/DemoRobertson/demorobertson.pdf).
Using the `@ode_def` macro from ParameterizedFunctions.jl, we can do this by writing:

```julia
using ParameterizedFunctions
f = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃
```

This looks just like pseudocode! The macro will expand this to the "standard form",
i.e. the ugly computer form:

```julia
function f(du,u,p,t)
  du[1] = -p[1]*u[1] + p[3]*u[2]*u[3]
  du[2] = p[1]*u[1] - p[2]*u[2]^2 - p[3]*u[2]*u[3]
  du[3] = p[2]*u[2]^2
end
```

Note that one doesn't need to use numbered variables: DifferentialEquations.jl
will number the variables for you. For example, the following defines the function
for the Lotka-Volterra model, with full Unicode support to boot:

```julia
f = @ode_def begin
  d🐁  = α*🐁  - β*🐁*🐈
  d🐈 = -γ*🐈 + δ*🐁*🐈
end α β γ δ
```

### Limitations

The macro is a Domain-Specific Language (DSL) and thus has different internal
semantics than standard Julia functions. In particular:

1) Control sequences and conditionals (while, for, if) will not work in the macro.
2) Intermediate calculations (likes that don't start with `d_`) are incompatible
   with the Jacobian etc. calculations.
3) The macro has to use `t` for the independent variable.

### Extra Optimizations

Because the ParameterizedFunction defined by the macro holds the definition at a
symbolic level, optimizations are provided by SymEngine. Using the symbolic
calculator, in-place functions for many things such as Jacobians, Hessians, etc.
are symbolically pre-computed. In addition, functions for the inverse Jacobian,
Hessian, etc. are also pre-computed. In addition, parameter gradients and
Jacobians are also used.

Normally these will be computed fast enough that the user doesn't have to worry.
However, in some cases you may want to restrict the number of functions (or get rid
of a warning). For more information,
please see the [ParameterizedFunctions.jl documentation](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl).
