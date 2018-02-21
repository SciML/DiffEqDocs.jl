# Immutable `struct` Handling

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
modified_struct = remake(original_struct;
  field_1 = value_1,
  field_2 = value_2,
  ...
)
```

where `field_N` and `value_N` are renamed to appropriate field names
and new desired values.
