# Dynamical, Hamiltonian and 2nd Order ODE Problems

## Mathematical Specification of a Dynamical ODE Problem

These algorithms require a Partitioned ODE of the form:

```math
\frac{du}{dt} = f_1(v) \\
\frac{dv}{dt} = f_2(u,p,t) \\
```
This is a Partitioned ODE partitioned into two groups, so the functions should be
specified as `f1(du,u,v,p,t)` and `f2(dv,u,v,p,t)` (in the inplace form), where `f1`
is independent of `t` and `u`, and unless specified by the solver,
`f2` is independent of `v`. This includes discretizations arising from
`SecondOrderODEProblem`s where the velocity is not used in the acceleration function,
and Hamiltonians where the potential is (or can be) time-dependent but the kinetic
energy is only dependent on `v`.

Note that some methods assume that the integral of `f1` is a quadratic form. That
means that `f1=v'*M*v`, i.e. ``\int f_1 = \frac{1}{2} m v^2``, giving `du = v`.
This is equivalent to saying that the kinetic energy is related to ``v^2``. The
methods which require this assumption will lose accuracy if this assumption is
violated. Methods listed make note of this requirement with "Requires
quadratic kinetic energy".

### Constructor

```julia
DynamicalODEProblem{isinplace}(f1,f2,u0,v0,tspan,callback=CallbackSet(),mass_matrix=I)
```

Defines the ODE with the specified functions. `isinplace` optionally sets whether
the function is inplace or not. This is determined automatically, but not inferred.

### Fields

* `f1` and `f2`: The functions in the ODE.
* `u0`: The initial condition.
* `du0`: The initial derivative.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

## Mathematical Specification of a 2nd Order ODE Problem

To define a 2nd Order ODE Problem, you simply need to give the function ``f``
and the initial condition ``u₀`` which define an ODE:

```math
u'' = f(u,u',p,t)
```

`f` should be specified as `f(u,du,p,t)` (or in-place as `f(ddu,u,du,p,t)`), and `u₀`
should be an AbstractArray (or number) whose geometry matches the desired
geometry of `u`. Note that we are not limited to numbers or vectors for `u₀`;
one is allowed to provide `u₀` as arbitrary matrices / higher dimension tensors
as well.

From this form, a dynamical ODE:

```math
u' = v \\
v' = f(u,p,t,v) \\
```

is generated.

### Constructors

```julia
SecondOrderODEProblem{isinplace}(f,u0,du0,tspan,callback=CallbackSet(),mass_matrix=I)
```

Defines the ODE with the specified functions.

### Fields

* `f`: The function in the ODE.
* `u0`: The initial condition.
* `du0`: The initial derivative.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.

## Hamiltonian Problems

`HamiltonianProblem`s are provided by DiffEqPhysics.jl and provide an easy way
to define equations of motion from the corresponding Hamiltonian. To define a
`HamiltonianProblem` one only needs to specify the Hamiltonian:

```math
H(q,p)
```

and autodifferentiation (via ForwardDiff.jl) will create the appropriate
equations.

### Constructors

```julia
HamiltonianProblem{T}(H,q0,p0,tspan;kwargs...)
```

### Fields

* `H`: The Hamiltonian `H(q,p,params)` which returns a scalar.
* `q0`: The initial positions.
* `p0`: The initial momentums.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem.
  Defaults to nothing.
* `mass_matrix`: The mass-matrix. Defaults to `I`, the `UniformScaling` identity matrix.
