# ODE Problems

## Mathematical Specification of an ODE Problem

To define an ODE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
\frac{du}{dt} = f(t,u)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

`ODEProblem(f,u0,tspan,callback=CallbackSet(),mass_matrix=I)` : Defines the ODE with the specified functions.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

# Refined ODE Problems

The refined ODE types are types that specify the ODE to a much greater degree of
detail, and thus give the solver more information and make it easier to optimize.
There are three different kinds of refined problems: split (IMEX) problems,
partitioned problems, and constrained problems.

## Mathematical Specification of a Split ODE Problem

To define a `ODEProblem` in split form, you simply need to give a tuple of
functions ``(f_1,f_2,\ldots,f_n)`` and the initial condition ``u₀`` which
define an ODE:

```math
\frac{du}{dt} =  f_1(t,u) + f_2(t,u) + \ldots + f_n(t,u)
```

`f` should be specified as `f(t,u)` (or in-place as `f(t,u,du)`), and `u₀` should
be an AbstractArray (or number) whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Mathematical Specification of a Partitioned ODE Problem

To define a Partitioned `ODEProblem`, you need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the tuple of initial conditions ``(u₀,v₀,...)``
(tuple of the same size) which define an ODE:

```math
\frac{du}{dt} = f_1(t,u,v,...) \\
\frac{dv}{dt} = f_2(t,u,v,...) \\
...
```
like e.g. `ODEProblem((f_1, f_2,...), (u0, v0, ...), tspan)`. Each of the `f_i`
should be specified as `f_1(t,u,v,...)`, or in-place as `f_1(t,u,v,...,du)`.
The initial conditions should be `AbstractArray`s (or numbers) whose geometry matches
the desired geometry of `u`. Note that we are not limited to numbers or vectors
for `u₀`; one is allowed to provide `u₀` as arbitrary matrices / higher dimension
tensors as well. 

In some cases, the solvers may specify the functions in a split
form, for example:

```math
\frac{du}{dt} = f_1(t,u,v,...) + f_2(t,u,v,...) \\
\frac{dv}{dt} = f_3(t,u,v,...) \\
...
```
see the [solver's documentation](solvers/ode_solve) for the form it is expecting.

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

## Example Problems

Example problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/ode_premade_problems.jl).

To use a sample problem, such as `prob_ode_linear`, you can do something like:

```julia
# Pkg.add("DiffEqProblemLibrary")
using DiffEqProblemLibrary
prob = prob_ode_linear
sol = solve(prob)
```

```@docs
DiffEqProblemLibrary.prob_ode_linear
DiffEqProblemLibrary.prob_ode_2Dlinear
DiffEqProblemLibrary.prob_ode_bigfloatlinear
DiffEqProblemLibrary.prob_ode_bigfloat2Dlinear
DiffEqProblemLibrary.prob_ode_large2Dlinear
DiffEqProblemLibrary.prob_ode_2Dlinear_notinplace
DiffEqProblemLibrary.prob_ode_threebody
DiffEqProblemLibrary.prob_ode_pleides
DiffEqProblemLibrary.prob_ode_vanderpol
DiffEqProblemLibrary.prob_ode_vanderpol_stiff
DiffEqProblemLibrary.prob_ode_rober
DiffEqProblemLibrary.prob_ode_rigidbody
```
