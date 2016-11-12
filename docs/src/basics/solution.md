# Solution Handling

## Accessing the Values

The solution type has a lot of built in functionality to help analysis. For example,
it has an array interface for accessing the values. Internally, the solution type
has two important fields:

1. `u` which holds the Vector of values at each timestep
2. `t` which holds the times of each timestep.

Different solution types may add extra information as necessary, such as the
derivative at each timestep `du` or the spatial discretization `x`, `y`, etc.

Instead of working on the `Vector{uType}` directly, we can use the provided
array interface.

```julia
sol[i]
```

to access the value at timestep `i` (if the timeseres was saved), and

```julia
sol.t[i]
```

to access the value of `t` at timestep `i`. For multi-dimensional systems, this
will address first by time and secondly by component, and thus

```julia
sol[i,j]
```

will be the `j`th component at timestep `i`. If the independent variables had shape
(for example, was a matrix), then `j` is the linear index. We can also access
solutions with shape:

```julia
sol[i,j,k]
```

gives the `[j,k]` component of the system at timestep `i`. The colon operator is
supported, meaning that

```julia
sol[:,j]
```

gives the timeseries for the `j`th component.

If the solver allows for dense output and `dense=true` was set for the solving
(which is the default), then we can access the approximate value
at a time `t` using the command

```julia
sol(t)
```

Note that the interpolating function allows for `t` to be a vector and use this
to speedup the interpolation calculations.

The solver interface also gives tools for using comprehensions over the solution.
Using the `tuples(sol)` function, we can get a tuple for the output at each
timestep. This allows one to do the following:

```julia
[t+2u for (t,u) in tuples(sol)]
```

One can use the extra components of the solution object as well using `zip`. For
example, say the solution type holds `du`, the derivative at each timestep. One
can comprehend over the values using:

```julia
[t+3u-du for (t,u,du) in zip(sol.t,sol.u,sol.du)]
```

Note that the solution object acts as a vector in time, and so its length is the
number of saved timepoints.

## Special Fields

The solution interface also includes some special fields. The problem object
`prob` and the algorithm used to solve the problem `alg` are included in the
solution. Additionally, the field `dense` is a boolean which states whether
the interpolation functionality is available. Lastly, there is a mutable state
`tslocation` which controls the plot recipe behavior. By default, `tslocation=0`.
Its values have different meanings between partial and ordinary differntial equations:

- `tslocation=0`  for non-spatial problems (ODEs) means that the plot recipe
  will plot the full solution. `tslocation=i` means that it will only plot the
  timepoint `i`.
- `tslocation=0` for spatial problems (PDEs) means the plot recipe will plot
  the final timepoint. `tslocation=i` means that the plot recipe will plot the
  `i`th timepoint.

What this means is that for ODEs, the plots will default to the full plot and PDEs
will default to plotting the surface at the final timepoint. The iterator interface
simply iterates the value of `tslocation`, and the `animate` function iterates
the solution calling solve at each step.

## Problem-Specific Features

Extra fields for solutions of specific problems are specified in the appropriate
problem definition page.  
