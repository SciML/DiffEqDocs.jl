# Refined DAE Problems

The refined DAE types are types that specify the DAE to a much greater degree of
detail, and thus give the solver more information and make it easier to optimize.
There are three different kinds of refined problems: split (IMEX) problems,
partitioned problems, and constrained problems.

## Mathematical Specification of a Split DAE Problem

To define a `SplitDAEProblem`, you simply need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the initial condition ``u₀`` which define an ODE:

```math
0 = f_1(t,u,u') + f_2(t,u,u') + \ldots + f_n(t,u,u')
```

`f` should be specified as `f(t,u,du)` (or in-place as `f(t,u,du,res)`), and `u₀`
should be an AbstractArray (or number) whose geometry matches the desired geometry
of `u`.

### Constructors

`SplitDAEProblem(f,u0,tspan,callback=nothing,mass_matrix=I)` : Defines the ODE with the specified functions.

### Fields

* `f`: The tuple of functions in the ODE.
* `u0`: The initial condition.
* `du0`: The initial derivative condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

### Special Solver Options

### Special Solution Fields

None. It returns a standard DAE solution.

## Mathematical Specification of a Partitioned ODE Problem

To define a `PartitionedDAEProblem`, you need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the tuple of initial conditions ``(u₀,v₀,...)`` (tuple
of the same size) which define an ODE:

```math
\frac{du}{dt} = f_1(t,u,v,...,du,dv,...) \\
\frac{dv}{dt} = f_2(t,u,v,...,du,dv,...) \\
...
```

`f` should be specified as `f(t,u,v,...,du,dv,...)` (or in-place as `f(t,u,v,...,du,dv,...,res)`), and
the initial conditions should be AbstractArrays (or numbers) whose geometry matches
the desired geometry of `u`. Note that we are not limited to numbers or vectors
for `u₀`; one is allowed to provide `u₀` as arbitrary matrices / higher dimension
tensors as well.

### Constructors

`PartitionedDAEProblem(f,u0,tspan,callback=nothing,mass_matrix=I)` : Defines the ODE with
the specified functions.

### Fields

* `f`: The tuple of functions for the ODE.
* `u0`: The tuple of initial conditions.
* `du0`: The tuple of initial derivatives.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

### Special Solver Options

### Special Solution Fields
