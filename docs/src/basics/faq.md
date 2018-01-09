# Frequently Asked Questions

This page is a compilation of frequently asked questions and answers.

## Complicated Models

#### Can I switch my ODE function in the middle of integration?

There are a few ways to do this. The simplest way is to just have a parameter to
switch between the two. For example:

```julia
function pf_func(t,u,p,du)
  if p == 0
    du[1] = 2u[1]
  else
    du[1] = -2u[1]
  end
  du[2] = -u[2]
end
pf = ParameterizedFunction(pf_func,0)
```

Then in a callback you can make the `affect!` function modify `integrator.prob.f.params`.
For example, we can make it change when `u[2]<0.5` via:

```julia
condition(t,u,integrator) = u[2] - 0.5
affect!(integrator) = integrator.prob.f.params = 1
```

Then it will change betweeen the two ODE choices for `du1` at that moment.
Another way to do this is to make the ODE functions all be the same type
via FunctionWrappers.jl, but that is unnecessary. With the way that modern
processors work, there exists branch prediction and thus execution of a conditional
is free if it's predictable which branch will be taken. In this case, almost every
call to `f` takes the `p==0` route until the callback, at which point it is
almost always the `else` route. Therefore the processor will effectively get
rid of the computational cost associated with this, so you're likely over-optimizing
if you're going further (unless this change happens every step, but even then
this is probably the cheapest part of the computation...).

## Numerical Error

#### The solver doesn't obey physical law X (e.g. conservation of energy)

Yes, this is because the numerical solution of the ODE is not the exact solution.
There are a few ways that you can handle this problem. One way is to get a more
exact solution. Thus instead of

```julia
sol = solve(prob,alg)
```

use

```julia
sol = solve(prob,alg,abstol=1e-10,reltol=1e-10)
```

Of course, there's always a tradeoff between accuracy and efficiency, so play
around to find out what's right for your problem.

Another thing you can do is use a callback. There are some
[premade callbacks in the callback library](http://docs.juliadiffeq.org/latest/features/callback_library.html) which
handle these sorts of things like projecting to manifolds and preserving positivity.

##### The symplectic integrator doesn't conserve energy?

Yes, symplectic integrators do not exactly conserve energy. It is a common
misconception that they do. What symplectic integrators actually do is solve
for a trajectory which rests on a symplectic manifold that is perturbed from
the true solution's manifold by the truncation error. This means that symplectic
integrators do not experience (very much) long time drift, but their orbit is
not exactly the same as the true solution in phase space and thus you will
see differences in energy that tend to look periodic. There is a small drift
which grows linearly and is related to floating point error, but this drift
is much less than standard methods. This is why symplectic methods are recommended
for long time integration.

For conserving energy, there are a few things you can do. First of all, the energy
error is related to the integration error, so simply solving with higher accuracy
will reduce the error. The results in the
[DiffEqBenchmarks](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl) show
that using a `DPRKN` method with low tolerance can be a great choice. Another
thing you can do is use
[the ManifoldProjection callback from the callback library](http://docs.juliadiffeq.org/latest/features/callback_library.html).

## Autodifferentiation and Dual Numbers

#### Are the native Julia solvers compatible with autodifferentiation?

Yes! If the algorithm does not use adaptive time stepping, then you simply
need to make the initial condition have elements of Dual numbers. If the
algorithm uses Dual numbers, you need to make sure that time is also
given by Dual numbers. A quick explanation of this is because changing
the value of the initial condition will change the error in the steps, thus
causing different steps to be taken changing the time values.

To show this in action, let's say we want to find the Jacobian of solution
of the Lotka-Volterra equation at `t=10` with respect to the parameters.

```julia
function pf_func(t,u,p,du)
  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
end
function f(p)
  pf = ParameterizedFunction(pf_func,p)
  prob = ODEProblem(pf,eltype(p).([1.0,1.0]),eltype(p).((0.0,10.0)))
  solve(prob,Tsit5(),save_everystep=false)[end]
end
```

This function takes in new parameters and spits out the solution at the end.
We make the inital condition `eltype(p).([1.0,1.0])` so that way it's typed to
be Dual numbers whenever `p` is an array of `Dual` numbers, and we do the same
for the timespan. Then we can take the Jacobian via ForwardDiff.jl:

```julia
using ForwardDiff
ForwardDiff.jacobian(f,[1.5,1.0])

2×2 Array{Float64,2}:
  2.20214   0.189782
 -6.2273   -0.700188
```

and compare it to Calculus.jl:

```julia
Calculus.jacobian(f,[1.5,1.0],:central)

2×2 Array{Float64,2}:
  2.20214   0.189782
 -6.2273   -0.700188
```

#### I get Dual number errors when I solve my ODE with Rosenbrock or SDIRK methods...?

This is because you're using a cache which is not compatible with autodifferentiaion
via ForwardDiff.jl. For example, if we use the ODE function:

```julia
const tmp = zeros(4)
const A = rand(4,4)
function f(t,u,du)
  A_mul_B!(tmp,A,u)
  du .= tmp .+ u
end
```

Here we use a cached temporary array in order to avoid the allocations of matrix
multiplication. When autodifferentiation occurs, the element type of `u` is
`Dual` numbers, so `A*u` produces `Dual` numbers, so the error arises when it
tries to write into `tmp`. There are two ways to avoid this. The first way,
the easy way, is to just turn off autodifferentiation with the `autodiff=false`
option in the solver. Every solver which uses autodifferentiation has this option.
Thus we'd solve this with:

```julia
prob = ODEProblem(f,rand(4),(0.0,1.0))
sol = solve(prob,Rosenbrock23(autodiff=false))
```

and it will use a numerical differentiation fallback (DiffEqDiffTools.jl) to
calculate Jacobians.

** Warning: Advanced **

The more difficult way is to create a Dual cache which dispatches for the cache
choice based on the element type of `u`. This is done by the following:

```julia
using ForwardDiff
struct MyTag end
immutable DiffCache{T<:AbstractArray, S<:AbstractArray}
    du::T
    dual_du::S
end

function DiffCache{chunk_size}(T, size, ::Type{Val{chunk_size}})
    DiffCache(zeros(T, size...), zeros(ForwardDiff.Dual{nothing,T,chunk_size}, size...))
end

DiffCache(u::AbstractArray) = DiffCache(eltype(u),size(u),Val{ForwardDiff.pickchunksize(length(u))})

get_tmp{T<:ForwardDiff.Dual}(dc::DiffCache, ::Type{T}) = dc.dual_du
get_tmp(dc::DiffCache, T) = dc.du
```

Now we can get a cache that by dispatch either gives a cache array of `Dual`
numbers or just floating point numbers:


```julia
const dual_cache = DiffCache(rand(4)) # Build the cache, this must match your IC
du = get_tmp(dual_cache,typeof(rand(4))) # Gives a Array{Float64}
dual_du = get_tmp(dual_cache,typeof(ForwardDiff.Dual(0.2,3.0))) # Gives Array{Dual}
```

Note that you have to make sure that your chunk size matches the choice in the
ODE solver (by default, it uses `ForwardDiff.pickchunksize(length(u))` as well,
so you only need to change this if you explicitly set `chunksize = ...`). Now
we can setup and solve our ODE using this cache:

```julia
function f(t,u,du)
  # Get du from cache
  tmp = get_tmp(dual_cache,eltype(u))
  # Fix tag
  _tmp = reinterpret(eltype(u),tmp)
  A_mul_B!(_tmp,A,u)
  du .= _tmp .+ u
end
prob = ODEProblem(f,rand(4),(0.0,1.0))
sol = solve(prob,Rosenbrock23())
```

Small explanation is in order. `tmp = get_tmp(dual_cache,eltype(u))` makes `tmp`
match `u` in terms of `Dual` or not, but it doesn't necessarily match the tag.
So now we `reinterpret` our `Dual` array to put the right tag on there. Note
that this simply changes type information and thus does not create any temporary
arrays. Once we do that, our cached array is now typed correctly to hold the
result of `A_mul_B!`.
