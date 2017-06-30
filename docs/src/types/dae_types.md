# DAE Problems

## Mathematical Specification of an DAE Problem

To define a DAE Problem, you simply need to give the function ``f`` and the initial
condition ``u₀`` which define an ODE:

```math
0 = f(t,u,du)
```

`f` should be specified as `f(t,u,du)` (or in-place as `f(t,u,du,resid)`).
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well.

## Problem Type

### Constructors

`DAEProblem(f,u0,du0,tspan)` : Defines the ODE with the specified functions.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `du0`: The initial condition for the derivative.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to a black CallbackSet, which will have no effect.
* `differential_vars`: A logical array which declares which variables are the
  differential (non algebraic) vars (i.e. `du'` is in the equations for this
  variable). Defaults to nothing. Some solvers may require this be set if an
  initial condition needs to be determined.

# Refined DAE Problems

The refined DAE types are types that specify the DAE to a much greater degree of
detail, and thus give the solver more information and make it easier to optimize.
There are three different kinds of refined problems: split (IMEX) problems,
partitioned problems, and constrained problems.

## Mathematical Specification of a Split DAE Problem

To define a split `DAEProblem`, you simply need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the initial condition ``u₀`` which define an ODE:

```math
0 = f_1(t,u,u') + f_2(t,u,u') + \ldots + f_n(t,u,u')
```

`f` should be specified as `f(t,u,du)` (or in-place as `f(t,u,du,res)`), and `u₀`
should be an AbstractArray (or number) whose geometry matches the desired geometry
of `u`.

## Mathematical Specification of a Partitioned ODE Problem

To define a `PartitionedDAEProblem`, you need to give a tuple of functions
``(f_1,f_2,\ldots,f_n)`` and the tuple of initial conditions ``(u₀,v₀,...)``
(tuple of the same size) which define an ODE:

```math
\frac{du}{dt} = f_1(t,u,v,...,du,dv,...) \\
\frac{dv}{dt} = f_2(t,u,v,...,du,dv,...) \\
```

`f` should be specified as `f(t,u,v,...,du,dv,...)` (or in-place as
`f(t,u,v,...,du,dv,...,res)`), and the initial conditions should be
AbstractArrays (or numbers) whose  geometry matches
the desired geometry of `u`. Note that we are not limited to numbers or vectors
for `u₀`; one is allowed to provide `u₀` as arbitrary matrices / higher dimension
tensors as well.

## Example Problems

Examples problems can be found in [DiffEqProblemLibrary.jl](https://github.com/JuliaDiffEq/DiffEqProblemLibrary.jl/blob/master/src/dae_premade_problems.jl).

To use a sample problem, such as `prob_dae_resrob`, you can do something like:

```julia
#Pkg.add("DiffEqProblemLibrary")
using DiffEqProblemLibrary
prob = prob_dae_resrob
sol = solve(prob,IDA())
```
