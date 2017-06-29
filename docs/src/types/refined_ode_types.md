# Refined ODE Problems

The refined ODE types are types that specify the ODE to a much greater degree of
detail, and thus give the solver more information and make it easier to optimize.
There are three different kinds of refined problems: split (IMEX) problems,
partitioned problems, and constrained problems.

## Mathematical Specification of a Split ODE Problem

To define a `SplitODEProblem`, you simply need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the initial condition ``u₀`` which define an ODE:

```math
\frac{du}{dt} =  f_1(t,u) + f_2(t,u) + \ldots + f_n(t,u)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

### Constructors

`SplitODEProblem(f,u0,tspan,callback=nothing,mass_matrix=I)` : Defines the ODE with the specified functions.

### Fields

* `f`: The tuple of functions in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

## Mathematical Specification of a Partitioned ODE Problem

To define a `PartitionedODEProblem`, you need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the tuple of initial conditions ``(u₀,v₀,...)`` (tuple
of the same size) which define an ODE:

```math
\frac{du}{dt} = f_1(t,u,v,...) \\
\frac{dv}{dt} = f_2(t,u,v,...) \\
...
```

`f` should be specified as `f(t,u,v,...)` (or in-place as `f(t,u,v,...,du)`), and
the initial conditions should be AbstractArrays (or numbers) whose geometry matches
the desired geometry of `u`. Note that we are not limited to numbers or vectors
for `u₀`; one is allowed to provide `u₀` as arbitrary matrices / higher dimension
tensors as well. In some cases, the solvers may specify the functions in a split
form, for example:

```math
\frac{du}{dt} = f_1(t,u,v,...) + f_2(t,u,v,...) \\
\frac{dv}{dt} = f_3(t,u,v,...) \\
...
```

See the solver's documentation for the form it is expecting.

### Constructors

`PartitionedODEProblem(f,u0,tspan,callback=nothing,mass_matrix=I)` : Defines the ODE with
the specified functions.

### Fields

* `f`: The tuple of functions for the ODE.
* `u0`: The tuple of initial conditions.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

## Mathematical Specification of an Second Order ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
u'' = f(t,u,u')
```

`f` should be specified as `f(t,u,du)` (or in-place as `f(t,u,du,ddu)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

From this form, a partitioned ODE

```math
u' = v \\
v' = f(t,u,v) \\
```

is generated.

### Problem Type

### Constructors

`SecondOrderODEProblem(f,u0,du0,tspan,callback=CallbackSet(),mass_matrix=I)` : Defines the ODE with the specified functions.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `du0`: The initial derivative.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.
