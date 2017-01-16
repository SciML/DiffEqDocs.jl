# ParameterizedFunctions

## Transforming User-Defined Functions to ParameterizedFunctions

### ParameterizedFunction Constructor

The easiest way to make a `ParameterizedFunction` is to use the constructor:

```julia
pf = ParameterizedFunction(f,params)
```

The form for `f` is `f(t,u,params,du)`
where `params` is any type which defines the parameters. The
resulting `ParameterizedFunction` has the function call `pf(t,u,params,du)`
which matches the original function, and a call `pf(t,u,du)` which uses internal
parameters which can be used with a differential equation solver. Note that the
internal parameters can be modified at any time via the field: `pf.p = ...`.

An additional version exists for `f(t,u,params)` which will then act as the
not in-place version `f(t,u)` in the differential equation solvers.

Note that versions exist for the other types of differential equations as well.
There are

```julia
pf = DAEParameterizedFunction(f,params)
pf = DDEParameterizedFunction(f,params)
```

for DAEs and DDEs respectively. For DAEs, the in-place syntax is `f(t,u,params,du,out)`
and the not in-place syntax is `f(t,u,params,du)`. For DDEs, the in-place syntax is
`f(t,u,h,params,du)` and the not in-place syntax is `f(t,u,h,params)`

### Examples

```julia
pf_func = function (t,u,p,du)
  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
end

pf = ParameterizedFunction(pf_func,[1.5,1.0])
```

And now `pf` can be used in the differential equation solvers and the ecosystem
functionality which requires explicit parameters (parameter estimation, etc.).

Note that the not in-place version works the same:

```julia
pf_func2 = function (t,u,p)
  [p[1] * u[1] - p[2] * u[1]*u[2];-3 * u[2] + u[1]*u[2]]
end

pf2 = ParameterizedFunction(pf_func2,[1.5,1.0])
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
f = @ode_def ROBERExample begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁=>0.04 k₂=>3e7 k₃=>1e4
```

This looks just like pseudocode! The macro will expand this to the "standard form",
i.e. the ugly computer form:

```julia
f = (t,u,du) -> begin
  du[1] = -0.04*u[1] + 1e4*u[2]*u[3]
  du[2] = 0.04*u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3]
  du[3] = 3e7*u[2]^2
end
```

Note that one doesn't need to use numbered variables: DifferentialEquations.jl
will number the variables for you. For example, the following defines the function
for the Lotka-Volterra model:

```julia
f = @ode_def LotkaVolterraExample begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1.0 c=>3.0 d=1.0
```

## Extra Features

Functions defined using the `@ode_def` macro come with many other features. For
example, since we used `=>` for `a`, `b`, and `c`, these parameters are explicitly
saved. That is, one can do:

```julia
f.a = 0.2
```

to change the parameter `f` to `0.2`. We can create a new function with new parameters
using the name we gave the macro:

```julia
g = LotkaVolterraExample(a=0.3,b=20.3)
```

In this case, `c` will default to the value we gave it in the macro.

Since the parameters are explicit, these functions can be used to analyze how the
parameters affect the model. Thus ParameterizedFunctions, when coupled with the
solvers, forms the backbone of functionality such as parameter estimation, parameter
sensitivity analysis, and bifurcation analysis.

## Extra Optimizations

Because the ParameterizedFunction defined by the macro holds the definition at a
symbolic level, optimizations are provided by SymEngine. Using the symbolic
calculator, in-place functions for many things such as Jacobians, Hessians, etc.
are symbolically pre-computed. In addition, functions for the inverse Jacobian,
Hessian, etc. are also pre-computed. In addition, parameter gradients and
Jacobians are also used.

Normally these will be computed fast enough that the user doesn't have to worry.
However, in some cases you may want to restrict the number of functions (or get rid
of a warning). Macros like `@ode_def_nohes` turn off the Hessian calculations,
and `@ode_def_noinvjac` turns off the Jacobian inversion. For more information,
please see the [ParameterizedFunctions.jl documentation](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl).

## Finite Element Method Macros

The other macro which is currently provided is the `@fem_def` macro. This macro
is for parsing and writing FEM functions. For example, in the FEM methods you have
to use `x[:,1]` instead of `x` and `x[:,2]` instead of `y`. The macro will automatically
do this replacement, along with adding in parameters. Since FEM functions are more
general, we also have to give it the function signature. Using the macro looks like this:

```julia
f  = @fem_def (x) DataFunction begin
  sin(α.*x).*cos(α.*y)
end α=>π

a = 2π
b = 8π*π
gD = @fem_def (x) DirichletBC begin
  sin(α.*x).*cos(α.*y)/β
end α=>a β=>b
```

This is equivalent to the definition:

```julia
f(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])
gD(x) = sin(2π.*x[:,1]).*cos(2π.*x[:,2])/(8π*π)
```

The true power comes in when dealing with nonlinear equations. The second argument,
which we skipped over as `()`, is for listing the variables you wish to define the
equation by. Mathematically you may be using `u`,`v`,`w`, etc., but for array-based
algorithms you need to use `u[:,1]`,`u[:,2]`,etc. To avoid obfuscated code, the
`@fem_def` macro does this conversion. For example:

```julia
l = @fem_def (t,x,u) begin
  du = ones(length(u))-α*u
  dv = ones(length(v))-v
end α=>0.5
```
says there are two equations, one for `u`: (`ones(length(u))-α*u`) and one for `v`:
`(ones(length(v))-v)`. This expands to the equation:

```julia
l = (t,x,u)  -> [ones(size(x,1))-.5u[:,1]   ones(size(x,1))-u[:,2]]
```

When you have 10+ variables, using `@fem_def` leads to code which is much
easier to read!
