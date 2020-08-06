# [Non-autonomous Linear ODE / Lie Group Problems](@id dynamical_prob)

Non-autonomous linear ODEs show up in a lot of scientific problems where
the differential equation lives on a manifold such as Lie Group. In these
situtations, specialized solvers can be utilized to enforce physical bounds
on the solution and enhance the solving.

## Mathematical Specification of a Non-autonomous Linear ODE

These algorithms require a Non-autonomous linear ODE of the form:

```math
u^\prime = A(u,p,t)u
```

Where ``A`` is an [AbstractDiffEqOperator](@ref DiffEqOperators) that is 
multiplied against ``u``. Many algorithms specialize on the form of ``A``, 
such as ``A`` being a constant or ``A`` being only time-dependent (``A(t)``). 

### Construction

Creating a non-autonomous linear ODE is the same as an `ODEProblem`, except `f`
is represented by an `AbstractDiffEqOperator` (note: this means that any standard
ODE solver can also be applied to problems written in this form). As an example:

```julia
function update_func(A,u,p,t)
    A[1,1] = cos(t)
    A[2,1] = sin(t)
    A[1,2] = -sin(t)
    A[2,2] = cos(t)
end
A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (10, 50.))
```

defines a quasi-linear ODE ``u^\prime = A(t)u`` where the components of ``A`` are
the given functions. Using that formulation, we can see that the general form is
``u^\prime = A(u,p,t)u``, for example:

```julia
function update_func(A,u,p,t)
    A[1,1] = 0
    A[2,1] = 1
    A[1,2] = -2*(1 - cos(u[2]) - u[2]*sin(u[2]))
    A[2,2] = 0
end
```

has a state-dependent linear operator. Note that many other `AbstractDiffEqOperator`s
can be used and `DiffEqArrayOperator` is just one version that represents `A` via
a matrix (other choices are matrix-free).

Note that if ``A`` is a constant, then it is sufficient to supply ``A`` directly without
an `update_func`.

### Note About Affine Equations

Note that the affine equation

```math
u^\prime = A(u,p,t)u + g(u,p,t)
```

can be written as a linear form by extending the size of the system by one to have a
constant term of 1, and extending `A` to have a new row containing the values of
`g(u,p,t)`. In this way, these types of equations can be handled by these specialized
integrators.