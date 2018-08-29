# Dynamical, Hamiltonian and 2nd Order ODE Problems

## Mathematical Specification of a Dynamical ODE Problem

These algorithms require a Partitioned ODE of the form:

```math
\frac{dv}{dt} = f_1(u,t) \\
\frac{du}{dt} = f_2(v) \\
```
This is a Partitioned ODE partitioned into two groups, so the functions should be
specified as `f1(dv,v,u,p,t)` and `f2(du,v,u,p,t)` (in the inplace form), where `f1`
is independent of `v` (unless specified by the solver), and `f2` is independent
of `u` and `t`. This includes discretizations arising from
`SecondOrderODEProblem`s where the velocity is not used in the acceleration function,
and Hamiltonians where the potential is (or can be) time-dependent but the kinetic
energy is only dependent on `v`.

Note that some methods assume that the integral of `f2` is a quadratic form. That
means that `f2=v'*M*v`, i.e. ``\int f_2 = \frac{1}{2} m v^2``, giving `du = v`.
This is equivalent to saying that the kinetic energy is related to ``v^2``. The
methods which require this assumption will lose accuracy if this assumption is
violated. Methods listed make note of this requirement with "Requires
quadratic kinetic energy".

### Constructor

```julia
DynamicalODEProblem(f::DynamicalODEFunction,v0,u0,tspan,callback=CallbackSet())
DynamicalODEProblem{isinplace}(f1,f2,v0,u0,tspan,callback=CallbackSet())
```

Defines the ODE with the specified functions. `isinplace` optionally sets whether
the function is inplace or not. This is determined automatically, but not inferred.

### Fields

* `f1` and `f2`: The functions in the ODE.
* `v0` and `u0`: The initial conditions.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.

## Mathematical Specification of a 2nd Order ODE Problem

To define a 2nd Order ODE Problem, you simply need to give the function ``f``
and the initial condition ``u₀`` which define an ODE:

```math
u'' = f(u',u,p,t)
```

`f` should be specified as `f(du,u,p,t)` (or in-place as `f(ddu,du,u,p,t)`), and `u₀`
should be an AbstractArray (or number) whose geometry matches the desired
geometry of `u`. Note that we are not limited to numbers or vectors for `u₀`;
one is allowed to provide `u₀` as arbitrary matrices / higher dimension tensors
as well.

From this form, a dynamical ODE:

```math
v' = f(v,u,p,t) \\
u' = v \\
```

is generated.

### Constructors

```julia
SecondOrderODEProblem{isinplace}(f,du0,u0,tspan,callback=CallbackSet())
```

Defines the ODE with the specified functions.

### Fields

* `f`: The function for the second derivative.
* `du0`: The initial derivative.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.

## Hamiltonian Problems

`HamiltonianProblem`s are provided by DiffEqPhysics.jl and provide an easy way
to define equations of motion from the corresponding Hamiltonian. To define a
`HamiltonianProblem` one only needs to specify the Hamiltonian:

```math
H(p,q)
```

and autodifferentiation (via ForwardDiff.jl) will create the appropriate
equations.

### Constructors

```julia
HamiltonianProblem{T}(H,p0,q0,tspan;kwargs...)
```

### Fields

* `H`: The Hamiltonian `H(p,q,params)` which returns a scalar.
* `p0`: The initial momentums.
* `q0`: The initial positions.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
