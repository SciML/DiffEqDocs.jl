simultaneously# [Local Sensitivity Analysis (Automatic Differentiation)](@id sensitivity)

Sensitivity analysis, or automatic differentiation of the solver, is provided
by the DiffEq suite. The model sensitivities are the derivatives of the
solution ``u(t)`` with respect to the parameters. Specifically, the local
sensitivity of the solution to a parameter is defined by how much the solution
would change by changes in the parameter, i.e. the sensitivity of the ith
independent variable to the jth parameter is
`` \frac{\partial u_i}{\partial p_{j}}``.

Sensitivity analysis serves two major purposes. On one hand, the sensitivities
are diagnostics of the model which are useful for understand how
it will change in accordance to changes in the parameters. But another use is
simply because in many cases these derivatives are useful. Sensitivity analysis
provides a cheap way to calculate the gradient of the solution which can be
used in parameter estimation and other optimization tasks.

There are three types of sensitivity analysis. Local forward sensitivity
analysis directly gives the gradient of the solution with respect to each
parameter along the time series. The computational cost scales like `N*M`,
where `N` is the number of states and `M` is the number of parameters. While
this gives all of the information, it can be expensive for models with large
numbers of parameters. Local adjoint sensitivity analysis solves directly for
the gradient of some functional of the solution, such as a cost function or
energy functional, in a manner that is cheaper when the number of parameters is
large. Global Sensitivity Analysis methods are meant to be used for exploring the
sensitivity over a larger domain without calculating derivatives and are covered
on a different page.

## Installation

This functionality does not come standard with DifferentialEquations.jl.
To use this functionality, you must install DiffEqSensitivty.jl:

```julia
]add DiffEqSensitivity
using DiffEqSensitivity
```

## High Level Interface: `concrete_solve`

The highest level interface is provided by the function `concrete_solve`:

```julia
concrete_solve(prob,alg,u0=prob.u0,p=prob.p,args...;sensealg=InterpolatingAdjoint(),kwargs...)
```

This function is fully compatible with automatic differentiation libraries
like [Zygote.jl](https://github.com/FluxML/Zygote.jl) and will automatically
replace any calculations of the solution's derivative with a fast method.
The return of `concrete_solve` is a concretized version of `solve`, i.e. no
interpolations are possible but it has the same array-like semantics of the
full solution object. The keyword argument `sensealg` controls the dispatch
to the algorithm used for the sensitivity calculation. The choices are as follows:

- `InterpolatingAdjoint`: The default.

## Sensitivity Algorithms

## Efficiency of the Different Methods

For an analysis of which methods will be most efficient for computing the
solution derivatives for a given problem, consult our analysis
[in this arxiv paper](https://arxiv.org/abs/1812.01892). A general rule of thumb
is:

- Discrete Forward Sensitivity Analysis via ForwardDiff.jl is the fastest for
  ODEs with small numbers of parameters (<100).
- Adjoint senstivity analysis is the fastest when the number of parameters is
  sufficiently large. There are three configurations of note. Using
  `backsolve` is the fastest and uses the least memory, but is not
  guaranteed to be stable. Checkpointing is the slowest but uses O(1)
  memory and is stable. Interpolating is the second fastest, is stable,
  but requires the ability to hold the full forward solution and its
  interpolation in memory.
- The methods which use automatic differentiation support the full range of
  DifferentialEquations.jl features (SDEs, DDEs, events, etc.), but only work
  on native Julia solvers. The methods which utilize altered ODE systems only
  work on ODEs (without events), but work on any ODE solver.

#### Options

Options for handling the adjoint computation are set by passing a `SensitivityAlg`
type, e.g. `SensitivityAlg(backsolve=true)`. Additionally, if Gauss-Kronrod quadrature
is used, the options `ireltol` and `iabstol` into `adjoint_sensitivities` controls
the behavior of the quadrature. Example calls:

```julia
res = adjoint_sensitivities(sol,Rodas4(),dg,t,ireltol=1e-8)

res = adjoint_sensitivities(sol,Vern9(),dg,t,reltol=1e-8,
                            sensealg=SensitivityAlg(backsolve=true))
```

* `checkpointing`: When enabled, the adjoint solutions compute the Jacobians
  by starting from the nearest saved value in `sol` and computing forward.
  By default, this is false if `sol.dense==true`, i.e. if `sol` has
  its higher order interpolation then this is by default disabled.
* `quad`: Use Gauss-Kronrod quadrature to integrate the adjoint sensitivity
  integral. Disabling it can decrease memory usage but increase computation
  time, since post-solution quadrature can be more accurate with less points
  using the continuous function. Default is `true`.
* `backsolve`: Solve the differential equation backward to get the past values.
  Note that for chaotic or non-reversible systems, such as though that solve to
  a steady state, enabling this option can lead to wildly incorrect results.
  Enabling it can decrease memory usage but increase computation time. When it
  is set to `true`, `quad` will be automatically set to `false`. Default is `false`.
* `autodiff`: Use automatic differentiation in the internal sensitivity algorithm
  computations. This is the only option that is passed, this will flip `autojacvec`
  to false, since that requires reverse-mode AD, and thus finite differencing for the
  full Jacobian will be used. Default is `true`.
* `chunk_size`: Chunk size for forward mode differentiation if full Jacobians are
  built (`autojacvec=false` and `autodiff=true`). Default is `0` for automatic
  choice of chunk size.
* `autojacvec`: Calculate the vector-Jacobian (adjoint sensitivity analysis) product
  via automatic differentiation with special seeding. Due to being a `vjp` function,
  this option requires using automatic differentiation, currently implemented with
  Tracker.jl. Default is `true` if `autodiff` is true, and otherwise is false. If
  `autojacvec=false`, then a full Jacobian has to be computed, and this will default
  to using a `f.jac` function provided by the user from the problem of the forward pass.
  Otherwise, if `autodiff=true` then it will use forward-mode AD for the Jacobian, otherwise
  it will fall back to using a numerical approximation to the Jacobian.

A way to understand these options at a higher level is as follows:

* For the choice of the overall sensitivity calculation algorithm, using interpolation is
  preferred if the `sol` is continuous. Otherwise it will use checkpointed adjoints by default
  if the user passes in a `sol` without a good interpolation. Using `backsolve` is unstable
  except in specific conditions, and thus is only used when chosen by the user.
* In any of these cases `quad=false` is the default, which enlarges the ODE system to calculate
  the integral simultaniously to the ODE solution. This reduces the memory cost, though in some
  cases solving the reverse problem with a continuous solution and then using QuadGK.jl to
  perform the quadrature can use less reverse-pass points and thus decrease the computation
  time, though this is rare when the number of parameters is large.
* Using automatic differentiation everywhere is the preferred default. This means the `vjp`
  will be performed using reverse-mode AD with Tracker.jl and no Jacobian will be formed.
  If `autodiff=false`, then `autojacvec=false` is set since it assumes that user function
  is not compatible with any automatic differentiation. In this case, if a user-defined
  Jacobian function exists, this will be used. Otherwise this means that the
  `vjp` will be computed by forming the Jacobian with finite differences and then doing
  the matrix-vector product. As an intermediate option, one can set `autodiff=true` with
  `autojacvec=false` to compute the Jacobian with forward-mode AD and then perform the
  vector-Jacobian product using that matrix.

## Lower Level Sensitivity Analysis Interfaces

While the high level interface is sufficient for interfacing with automatic
differentiation, for example allowing compatibility with neural network libraries,
in some cases one may want more control over the way the sensitivities are
calculated in order to squeeze out every ounce of optimization. If you're that
user, then this section of the docs is for you.

## Local Forward Sensitivity Analysis via ODELocalSensitivityProblem

Local forward sensitivity analysis gives a solution along with a timeseries of
the sensitivities. Thus if one wishes to have a derivative at every possible
time point, directly utilizing the `ODELocalSensitivityProblem` can be more
efficient.

### ODELocalSensitivityProblem Syntax

`ODELocalSensitivityProblem` is similar to an `ODEProblem`, but takes an
`AbstractSensitivityAlgorithm` that describes how to append the forward sensitivity
equation calculation to the time evolution to simultaneously compute the derivative
of the solution with respect to parameters.

```julia
ODEForwardSensitivityProblem(f::DiffEqBase.AbstractODEFunction,u0,
                             tspan,p=nothing,
                             sensealg::AbstractForwardSensitivityAlgorithm = ForwardSensitivity();
                             kwargs...)
```

Once constructed, this problem can be used in `solve`. The solution can be
deconstructed into the ODE solution and sensitivities parts using the
`extract_local_sensitivities` function, with the following dispatches:

```julia
extract_local_sensitivities(sol, asmatrix::Val=Val(false)) # Decompose the entire time series
extract_local_sensitivities(sol, i::Integer, asmatrix::Val=Val(false)) # Decompose sol[i]
extract_local_sensitivities(sol, t::Union{Number,AbstractVector}, asmatrix::Val=Val(false)) # Decompose sol(t)
```

For information on the mathematics behind these calculations, consult
[the sensitivity math page](@ref sensitivity_math)

#### Example using an ODELocalSensitivityProblem

To define a sensitivity problem, simply use the `ODELocalSensitivityProblem` type
instead of an ODE type. For example, we generate an ODE with the sensitivity
equations attached for the Lotka-Volterra equations by:

```julia
function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + u[1]*u[2]
end

p = [1.5,1.0,3.0]
prob = ODELocalSensitivityProblem(f,[1.0;1.0],(0.0,10.0),p)
```

This generates a problem which the ODE solvers can solve:

```julia
sol = solve(prob,DP8())
```

Note that the solution is the standard ODE system and the sensitivity system combined.
We can use the following helper functions to extract the sensitivity information:

```julia
x,dp = extract_local_sensitivities(sol)
x,dp = extract_local_sensitivities(sol,i)
x,dp = extract_local_sensitivities(sol,t)
```

In each case, `x` is the ODE values and `dp` is the matrix of sensitivities
The first gives the full timeseries of values and `dp[i]` contains the time series of the
sensitivities of all components of the ODE with respect to `i`th parameter.
The second returns the `i`th time step, while the third
interpolates to calculate the sensitivities at time `t`. For example, if we do:

```julia
x,dp = extract_local_sensitivities(sol)
da = dp[1]
```

then `da` is the timeseries for ``\frac{\partial u(t)}{\partial p}``. We can
plot this

```julia
plot(sol.t,da',lw=3)
```

transposing so that the rows (the timeseries) is plotted.

![Local Sensitivity Solution](../assets/sensitivityplot.png)

Here we see that there is a periodicity to the sensitivity which matches
the periodicity of the Lotka-Volterra solutions. However, as time goes on the
sensitivity increases. This matches the analysis of Wilkins in Sensitivity
Analysis for Oscillating Dynamical Systems.

We can also quickly see that these values are equivalent to those given by
automatic differentiation and numerical differentiation through the ODE solver:

```julia
using ForwardDiff, Calculus
function test_f(p)
  prob = ODEProblem(f,eltype(p).([1.0,1.0]),eltype(p).((0.0,10.0)),p)
  solve(prob,Vern9(),abstol=1e-14,reltol=1e-14,save_everystep=false)[end]
end

p = [1.5,1.0,3.0]
fd_res = ForwardDiff.jacobian(test_f,p)
calc_res = Calculus.finite_difference_jacobian(test_f,p)
```

Here we just checked the derivative at the end point.

#### Internal representation

For completeness, we detail the internal representation. When using
ForwardDiffSensitivity, the representation is with `Dual` numbers under the
standard interpretation. The values for the ODE's solution at time `i` are the
`ForwardDiff.value.(sol[i])` portions, and the derivative with respect to
parameter `j` is given by `ForwardDiff.partials.(sol[i])[j]`.

When using ForwardSensitivity, the solution to the ODE are the first `n`
components of the solution. This means we can grab the matrix of solution
values like:

```julia
x = sol[1:sol.prob.indvars,:]
```

Since each sensitivity is a vector of derivatives for each function, the sensitivities
are each of size `sol.prob.indvars`. We can pull out the parameter sensitivities from
the solution as follows:

```julia
da = sol[sol.prob.indvars+1:sol.prob.indvars*2,:]
db = sol[sol.prob.indvars*2+1:sol.prob.indvars*3,:]
dc = sol[sol.prob.indvars*3+1:sol.prob.indvars*4,:]
```

This means that `da[1,i]` is the derivative of the `x(t)` by the parameter `a`
at time `sol.t[i]`. Note that all of the functionality available to ODE solutions
is available in this case, including interpolations and plot recipes (the recipes
will plot the expanded system).

## Adjoint Sensitivity Analysis via adjoint_sensitivities

Adjoint sensitivity analysis is used to find the gradient of the solution
with respect to some functional of the solution. In many cases this is used
in an optimization problem to return the gradient with respect to some cost
function. It is equivalent to "backpropagation" or reverse-mode automatic
differentiation of a differential equation.

Using `adjoint_sensitivities` directly let's you do three things. One it can
allow you to be more efficient, since the sensitivity calculation can be done
directly on a cost function, avoiding the overhead of building the derivative
of the full concretized solution. It can also allow you to be more efficient
by directly controlling the forward solve that is then reversed over. Lastly,
it allows one to define a continuous cost function on the continuous solution,
instead of just at discrete data points.

### Syntax

There are two forms. For discrete adjoints, the form is:

```julia
du0,dp = adjoint_sensitivities(sol,alg,dg,ts;sensealg=InterpolatingAdjoint(),
                               checkpoints=sol.t,kwargs...)
```

where `alg` is the ODE algorithm to solve the adjoint problem, `dg` is the jump
function, `sensealg` is the sensitivity algorithm, and `ts` is the time points
for data. `dg` is given by:

```julia
dg(out,u,p,t,i)
```

which is the in-place gradient of the cost functional `g` at time point `ts[i]`
with `u=u(t)`.

For continuous functionals, the form is:

```julia
du0,dp = adjoint_sensitivities(sol,alg,g,nothing,dg;sensealg=InterpolatingAdjoint(),
                               checkpoints=sol.t,,kwargs...)
```

for the cost functional

```julia
g(u,p,t)
```

with in-place gradient

```julia
dg(out,u,p,t)
```

If the gradient is omitted, i.e.

```julia
du0,dp = adjoint_sensitivities(sol,alg,g,nothing;kwargs...)
```

then it will be computed automatically using ForwardDiff or finite
differencing, depending on the `autodiff` setting in the `AbstractSensitivityAlgorithm`.
Note that the keyword arguments are passed to the internal ODE solver for
solving the adjoint problem.

### Example discrete adjoints on a cost function

In this example we will show solving for the adjoint sensitivities of a discrete
cost functional. First let's solve the ODE and get a high quality continuous
solution:

```julia
function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + u[1]*u[2]
end

p = [1.5,1.0,3.0]
prob = ODEProblem(f,[1.0;1.0],(0.0,10.0),p)
sol = solve(prob,Vern9(),abstol=1e-10,reltol=1e-10)
```

Now let's calculate the sensitivity of the L2 error against 1 at evenly spaced
points in time, that is:

```math
L(u,p,t)=\sum_{i=1}^{n}\frac{\Vert1-u(t_{i},p)\Vert^{2}}{2}
```

for ``t_i = 0.5i``. This is the assumption that the data is `data[i]=1.0`.
For this function, notice we have that:

```math
\begin{aligned}
dg_{1}&=1-u_{1} \\
dg_{2}&=1-u_{2}
\end{aligned}
```

and thus:

```julia
dg(out,u,i) = (out.=1.0.-u)
```

If we had data, we'd just replace `1.0` with `data[i]`. To get the adjoint
sensitivities, call:

```julia
res = adjoint_sensitivities(sol,Vern9(),dg,t,abstol=1e-14,
                            reltol=1e-14)
```

This is super high accuracy. As always, there's a tradeoff between accuracy
and computation time. We can check this almost exactly matches the
autodifferentiation and numerical differentiation results:

```julia
using ForwardDiff,Calculus
function G(p)
  tmp_prob = remake(prob,u0=convert.(eltype(p),prob.u0),p=p)
  sol = solve(tmp_prob,Vern9(),abstol=1e-14,reltol=1e-14,saveat=t)
  A = convert(Array,sol)
  sum(((1-A).^2)./2)
end
G([1.5,1.0,3.0])
res2 = ForwardDiff.gradient(G,[1.5,1.0,3.0])
res3 = Calculus.gradient(G,[1.5,1.0,3.0])
res4 = Flux.Tracker.gradient(G,[1.5,1.0,3.0])
res5 = ReverseDiff.gradient(G,[1.5,1.0,3.0])
```

and see this gives the same values.

### Example controlling adjoint method choices and checkpointing

In the previous examples, all calculations were done using the interpolating
method. This maximizes speed but at a cost of requiring a dense `sol`. If it
is not possible to hold a dense forward solution in memory, then one can use
checkpointing. For example:

```julia
ts = [0.0,0.2,0.5,0.7]
sol = solve(prob,Vern9(),saveat=ts)
```

Creates a non-dense solution with checkpoints at `[0.0,0.2,0.5,0.7]`. Now we
can do:

```julia
res = adjoint_sensitivities(sol,Vern9(),dg,t,
                            sensealg=InterpolatingAdjoint(checkpointing=true))
```

When grabbing a Jacobian value during the backwards solution, it will no longer
interpolate to get the value. Instead, it will start a forward solution at the
nearest checkpoint to build local interpolants in a way that conserves memory.
By default the checkpoints are at `sol.t`, but we can override this:

```julia
res = adjoint_sensitivities(sol,Vern9(),dg,t,
                            sensealg=InterpolatingAdjoint(checkpointing=true),
                            checkpoints = [0.0,0.5])
```

#### Applicability of Backsolve and Caution

When `backsolve` is applicable it is a fast method and requires the least memory.
However, one must be cautious because not all ODEs are stable under backwards integration
by the majority of ODE solvers. An example of such an equation is the Lorenz equation.
Notice that if one solves the Lorenz equation forward and then in reverse with any
adaptive time step and non-reversible integrator, then the backwards solution diverges
from the forward solution. As a quick demonstration:

```julia
using Sundials, DiffEqBase
function lorenz(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-12)
prob2 = ODEProblem(lorenz,sol[end],(100.0,0.0))
sol = solve(prob,Tsit5(),reltol=1e-12,abstol=1e-12)
@show sol[end]-u0 #[-3.22091, -1.49394, 21.3435]
```

Thus one should check the stability of the backsolve on their type of problem before
enabling this method. Additionally, using checkpointing with backsolve can be a
low memory way to stabilize it.

### Example continuous adjoints on an energy functional

In this case we'd like to calculate the adjoint sensitivity of the scalar energy
functional:

```math
G(u,p)=\int_{0}^{T}\frac{\sum_{i=1}^{n}u_{i}^{2}(t)}{2}dt
```

which is:

```julia
g(u,p,t) = (sum(u).^2) ./ 2
```

Notice that the gradient of this function with respect to the state `u` is:

```julia
function dg(out,u,p,t)
  out[1]= u[1] + u[2]
  out[2]= u[1] + u[2]
end
```

To get the adjoint sensitivities, we call:

```julia
res = adjoint_sensitivities(sol,Vern9(),g,nothing,dg,abstol=1e-8,
                                 reltol=1e-8,iabstol=1e-8,ireltol=1e-8)
```

Notice that we can check this against autodifferentiation and numerical
differentiation as follows:

```julia
function G(p)
  tmp_prob = remake(prob,p=p)
  sol = solve(tmp_prob,Vern9(),abstol=1e-14,reltol=1e-14)
  res,err = quadgk((t)-> (sum(sol(t)).^2)./2,0.0,10.0,abstol=1e-14,reltol=1e-10)
  res
end
res2 = ForwardDiff.gradient(G,[1.5,1.0,3.0])
res3 = Calculus.gradient(G,[1.5,1.0,3.0])
```
