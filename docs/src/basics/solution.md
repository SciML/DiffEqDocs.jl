# [Solution Handling](@id solution)

## Accessing the Values

The solution type has a lot of built in functionality to help analysis. For example,
it has an array interface for accessing the values. Internally, the solution type
has two important fields:

1. `u` which holds the Vector of values at each timestep
2. `t` which holds the times of each timestep.

Different solution types may add extra information as necessary, such as the
derivative at each timestep `du` or the spatial discretization `x`, `y`, etc.

## Array Interface

Instead of working on the `Vector{uType}` directly, we can use the provided
array interface.

```julia
sol[j]
```

to access the value at timestep `j` (if the timeseries was saved), and

```julia
sol.t[j]
```

to access the value of `t` at timestep `j`. For multi-dimensional systems, this
will address first by component and lastly by time, and thus

```julia
sol[i,j]
```

will be the `i`th component at timestep `j`. Hence, `sol[j][i] == sol[i, j]`. This is done because Julia is column-major, so the leading dimension should be contiguous in memory. If the independent variables had shape
(for example, was a matrix), then `i` is the linear index. We can also access
solutions with shape:

```julia
sol[i,k,j]
```

gives the `[i,k]` component of the system at timestep `j`. The colon operator is
supported, meaning that

```julia
sol[i,:]
```

gives the timeseries for the `i`th component.

## Using the AbstractArray Interface

The `AbstractArray` interface can be directly used. For example, for a vector
system of variables `sol[i,j]` is a matrix with rows being the variables and
columns being the timepoints. Operations like `sol'` will
transpose the solution type. Functionality written for `AbstractArray`s can
directly use this. For example, the Base `cov` function computes correlations
amongst columns, and thus:

```julia
cov(sol)
```

computes the correlation of the system state in time, whereas

```julia
cov(sol,2)
```

computes the correlation between the variables. Similarly, `mean(sol,2)` is the
mean of the variable in time, and `var(sol,2)` is the variance. Other statistical
functions and packages which work on `AbstractArray` types will work on the
solution type.

At anytime, a true `Array` can be created using `Array(sol)`.

## Interpolations and Calculating Derivatives

If the solver allows for dense output and `dense=true` was set for the solving
(which is the default), then we can access the approximate value
at a time `t` using the command

```julia
sol(t)
```

Note that the interpolating function allows for `t` to be a vector and uses this to speed up the interpolation calculations. The full API for the interpolations is

```julia
sol(t,deriv=Val{0};idxs=nothing,continuity=:left)
```

The optional argument `deriv` lets you choose the number `n` derivative to solve the interpolation for, defaulting with `n=0`. Note that most of the derivatives have not yet been implemented (though it's not hard, it just has to be done by hand for each algorithm. Open an issue if there's a specific one you need). `continuity` describes whether to satisfy left or right continuity when a discontinuity is saved. The default is `:left`, i.e. grab the value before the callback's change, but can be changed to `:right`. `idxs` allows you to choose the indices the interpolation should solve for. For example,

```julia
sol(t,idxs=1:2:5)
```

will return a `Vector` of length 3 which is the interpolated values at `t` for components `1`, `3`, and `5`. `idxs=nothing`, the default, means it will return every component. In addition, we can do

```julia
sol(t,idxs=1)
```

and it will return a `Number` for the interpolation of the single value. Note that this interpolation only computes the values which are requested, and thus it's much faster on large systems to use this rather than computing the full interpolation and using only a few values.

In addition, there is an inplace form:

```julia
sol(out,t,deriv=Val{0};idxs=nothing,continuity=:left)
```

which will write the output to `out`. This allows one to use pre-allocated vectors for the output to improve the speed even more.

## Comprehensions

The solver interface also gives tools for using comprehensions over the solution.
Using the `tuples(sol)` function, we can get a tuple for the output at each
timestep. This allows one to do the following:

```julia
[t+2u for (u,t) in tuples(sol)]
```

One can use the extra components of the solution object as well as using `zip`. For
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
the interpolation functionality is available. Further, the field `destats`
contains the internal statistics for the solution process such as the number
of linear solves and convergence failures. Lastly, there is a mutable state
`tslocation` which controls the plot recipe behavior. By default, `tslocation=0`.
Its values have different meanings between partial and ordinary differential equations:

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

## Differential Equation Solver Statistics (destats)

```@docs
DiffEqBase.DEStats
```

## [Return Codes (RetCodes)](@id retcodes)

The solution types have a `retcode` field which returns a symbol signifying the
error state of the solution. The retcodes are as follows:

- `:Default`: The solver did not set retcodes.
- `:Success`: The integration completed without erroring or the steady state solver
  from `SteadyStateDiffEq` found the steady state.
- `:Terminated`: The integration is terminated with `terminate!(integrator)`.
  Note that this may occur by using `TerminateSteadyState` from the callback
  library `DiffEqCallbacks`.
- `:MaxIters`: The integration exited early because it reached its maximum number
  of iterations.
- `:DtLessThanMin`: The timestep method chose a stepsize which is smaller than the
  allowed minimum timestep, and exited early.
- `:Unstable`: The solver detected that the solution was unstable and exited early.
- `:InitialFailure`: The DAE solver could not find consistent initial conditions.
- `:ConvergenceFailure`: The internal implicit solvers failed to converge.
- `:Failure`: General uncategorized failures or errors.

## Problem-Specific Features

Extra fields for solutions of specific problems are specified in the appropriate
problem definition page.  
