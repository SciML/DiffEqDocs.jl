# Problem Interface

This page defines the common problem interface. There are certain rules that
can be applied to any function definition, and this page defines those behaviors.

## In-place vs Out-of-Place Function Definition Forms

Every problem definition has an in-place and out-of-place form, commonly referred
throughout DiffEq as IIP (`isinplace`) and OOP (out of place). The in-place form
is a mutating form. For example, on ODEs, we have that `f!(du,u,p,t)` is the
in-place form which, as its output, mutates `du`. Whatever is returned is simply
ignored. Similarly, for OOP we have the form `du=f(u,p,t)` which uses the return.

Each of the problem types have that the first argument is the option mutating
argument. The DiffEqBase system will automatically determine the functional
form and place a specifier `isinplace` on the function to carry as type information
whether the function defined for this `DEProblem` is in-place. However, every
constructor allows for manually specifying the in-placeness of the function.
For example, this can be done at the problem level like:

```julia
ODEProblem{true}(f,u0,tspan,p)
```

which declares that `isinplace=true`. Similarly this can be done at the
DEFunction level. For example:

```julia
ODEFunction{true}(f,jac=myjac)
```

## Type Specifications

Throughout DifferentialEquations.jl, the types that are given in a problem are
the types used for the solution. If an initial value `u0` is needed for a problem,
then the state variable `u` will match the type of that `u0`. Similarly, if
time exists in a problem the type for `t` will be derived from the types of the
`tspan`. Parameters `p` can be any type and the type will be matching how it's
defined in the problem.

For internal matrices, such as Jacobians and Brownian caches, these also match
the type specified by the user. `jac_prototype` and `rand_prototype` can thus
be any Julia matrix type which is compatible with the operations that will be
performed.

## Functional and Condensed Problem Inputs

Note that the initial condition can be written as a function of parameters and
initial time:

```julia
u0(p,t0)
```

and be resolved before going to the solver. Additionally, the initial condition
can be a distribution from Distributions.jl, in which case a sample initial condition
will be taken each time `init` or `solve` is called.

In addition, `tspan` supports the following forms. The single value form `t`
is equivalent to `(zero(t),t)`. The functional form is allowed:

```julia
tspan(p)
```

which outputs a tuple.

### Examples

```julia
prob = ODEProblem((u,p,t)->u,(p,t0)->p[1],(p)->(0.0,p[2]),(2.0,1.0))
using Distributions
prob = ODEProblem((u,p,t)->u,(p,t)->Normal(p,1),(0.0,1.0),1.0)
```

## Lower Level `__init` and `__solve`

At the high level, known problematic problems will emit warnings before entering
the solver to better clarify the error to the user. The following cases are
checked if the solver is adaptive:

- Integer times warn
- Dual numbers must be in the initial conditions and timespans
- Measurements.jl values must be in the initial conditions and timespans

If there is an exception to these rules, please file an issue. If one wants to
go around the high level solve interface and its warnings, one can call `__init`
or `__solve` instead.

## Modification of problem types

Problem-related types in DifferentialEquations.jl are immutable.  This
helps, e.g., parallel solvers to efficiently handle problem types.

However, you may want to modify the problem after it is created.  For
example, to simulate it for longer timespan.  It can be done by the
`remake` function:

```julia
prob1 = ODEProblem((u,p,t) -> u/2, 1.0, (0.0,1.0))
prob2 = remake(prob1; tspan=(0.0,2.0))
```

A general syntax of `remake` is

```julia
modified_problem = remake(original_problem;
  field_1 = value_1,
  field_2 = value_2,
  ...
)
```

where `field_N` and `value_N` are renamed to appropriate field names
and new desired values.
