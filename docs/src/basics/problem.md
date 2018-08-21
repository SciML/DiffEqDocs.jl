# Problem interface

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
