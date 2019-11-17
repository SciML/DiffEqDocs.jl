# Local Sensitivity Analysis (Automatic Differentiation)

Sensitivity analysis for ODE models is provided by the DiffEq suite. The model
sensitivities are the derivatives of the solution with respect to the
parameters. Specifically, the local sensitivity of the solution to a parameter
is defined by how much the solution would change by changes in the parameter,
i.e. the sensitivity of the ith independent variable to the jth parameter is
`` \frac{\partial u}{\partial p_{j}}``.

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

## Efficiency of the Different Methods

For an analysis of which methods will be most efficient for computing the
solution derivatives for a given problem, consult our analysis
[in this arxiv paper](https://arxiv.org/abs/1812.01892). A general rule of thumb
is:

- Discrete Forward Sensitivity Analysis via ForwardDiff.jl is the fastest for
  ODEs with small numbers of parameters (<100)
- Adjoint senstivity analysis is the fastest when the number of parameters is
  sufficiently large. There are three configurations of note. Using
  `backsolve` is the fastest and uses the least memory, but is not
  guerenteed to be stable. Checkpointing is the slowest but uses O(1)
  memory and is stable. Interpolating is the second fastest, is stable,
  but requires the ability to hold the full forward solution and its
  interpolation in memory.
- The methods which use automatic differentiation support the full range of
  DifferentialEquations.jl features (SDEs, DDEs, events, etc.), but only work
  on native Julia solvers. The methods which utilize altered ODE systems only
  work on ODEs (without events), but work on any ODE solver.

## Local Forward Sensitivity Analysis

Local forward sensitivity analysis gives a solution along with a timeseries of
the sensitivities along the solution.

### Discrete Local Forward Sensitivity Analysis via ForwardDiff.jl

This method is the application of ForwardDiff.jl numbers to the ODE solver. This
is done simply by making the `u0` state vector a vector of Dual numbers, and
multiple dispatch then allows the internals of the solver to propagate the
derivatives along the solution.

#### Examples using ForwardDiff.jl

The easiest way to use ForwardDiff.jl for local forward sensitivity analysis
is to simply put the ODE solve inside of a function which you would like to
differentiate. For example, let's define the ODE system for the Lotka-Volterra
equations:

```julia
function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + u[1]*u[2]
end

p = [1.5,1.0,3.0,1.0]
u0 = [1.0;1.0]
prob = ODEProblem(f,u0,(0.0,10.0),p)
```

Let's say we wanted to get the derivative of the final value w.r.t. each of
the parameters. We can define the following function:

```julia
function test_f(p)
  _prob = remake(prob;u0=convert.(eltype(p),prob.u0),p=p)
  solve(_prob,Vern9(),save_everystep=false)[end]
end
```

What this function does is use the `remake` function
[from the Problem Interface page](http://docs.juliadiffeq.org/latest/basics/problem/#Modification-of-problem-types-1)
to generate a new ODE problem with the new parameters, solves it, and returns
the solution at the final time point. Notice that it takes care to make sure
that the type of `u0` matches the type of `p`. This is because ForwardDiff.jl
will want to use Dual numbers, and thus to propagate the Duals throughout
the solver we need to make sure the initial condition is also of the type of
Dual number. On this function we can call ForwardDiff.jl and it will return
the derivatives we wish to calculate:

```julia
using ForwardDiff
fd_res = ForwardDiff.jacobian(test_f,p)
```

If we would like to get the solution and the value at the time point at the
same time, we can use [DiffResults.jl](https://github.com/JuliaDiff/DiffResults.jl).
For example, the following uses a single ODE solution to calculate the value
at the end point and its parameter Jacobian:

```julia
using DiffResults
res = DiffResults.JacobianResult(u0,p) # Build the results object
DiffResults.jacobian!(res,p) # Populate it with the results
val = DiffResults.value(res) # This is the sol[end]
jac = DiffResults.jacobian(res) # This is dsol/dp
```

If we would like to get the time series, we can do so by seeding the dual
numbers directly. To do this, we use the `Dual` constructor. The first value
is the value of the parameter. The second is a tuple of the derivatives. For
each value we want to take the derivative by, we seed a derivative with a
1 in a unique index. For example, we can build our parameter vector like:

```julia
using ForwardDiff: Dual
struct MyTag end
p1dual = Dual{MyTag}(1.5, (1.0, 0.0, 0.0, 0.0))
p2dual = Dual{MyTag}(1.0, (0.0, 1.0, 0.0, 0.0))
p3dual = Dual{MyTag}(3.0, (0.0, 0.0, 1.0, 0.0))
p4dual = Dual{MyTag}(3.0, (0.0, 0.0, 0.0, 1.0))
pdual = [p1dual, p2dual, p3dual, p4dual]
```

or equivalently using the `seed_duals` convenience function:

```julia
function seed_duals(x::AbstractArray{V},::Type{T},
                    ::ForwardDiff.Chunk{N} = ForwardDiff.Chunk(x)) where {V,T,N}
  seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})
  duals = [Dual{T}(x[i],seeds[i]) for i in eachindex(x)]
end
pdual = seed_duals(p,MyTag)
```

Next we need to make our initial condition Dual numbers so that these propogate
through the solution. We can do this manually like:

```julia
u0dual = [Dual{MyTag}(1.0, (0.0, 0.0, 0.0, 0.0)),Dual{MyTag}(1.0, (0.0, 0.0, 0.0, 0.0))]
```

or use the same shorthand from before:

```julia
u0dual = convert.(eltype(pdual),u0)
```

Now we just use these Dual numbers to solve:

```julia
prob_dual = ODEProblem(f,u0,tspan,pdual)
sol_dual = solve(prob_dual,Tsit5(), saveat=0.2)
```

The solution is now in terms of Dual numbers. We can extract the derivatives
by looking at the partials of the duals in the solution. For example, `sol[1,end]`
is the Dual number for the `x` component at the end of the integration, and so
`sol[1,end].partial[i]` is `dx(t_end)/dp_i`.

### Local Forward Sensitivity Analysis via ODELocalSensitivityProblem

For this method local sensitivity is computed using the sensitivity ODE:

```math
\frac{d}{dt}\frac{\partial u}{\partial p_{j}}=\frac{\partial f}{\partial u}\frac{\partial u}{\partial p_{j}}+\frac{\partial f}{\partial p_{j}}=J\cdot S_{j}+F_{j}
```

where

```math
J=\left(\begin{array}{cccc}
\frac{\partial f_{1}}{\partial u_{1}} & \frac{\partial f_{1}}{\partial u_{2}} & \cdots & \frac{\partial f_{1}}{\partial u_{k}}\\
\frac{\partial f_{2}}{\partial u_{1}} & \frac{\partial f_{2}}{\partial u_{2}} & \cdots & \frac{\partial f_{2}}{\partial u_{k}}\\
\cdots & \cdots & \cdots & \cdots\\
\frac{\partial f_{k}}{\partial u_{1}} & \frac{\partial f_{k}}{\partial u_{2}} & \cdots & \frac{\partial f_{k}}{\partial u_{k}}
\end{array}\right)
```

is the Jacobian of the system,

```math
F_{j}=\left(\begin{array}{c}
\frac{\partial f_{1}}{\partial p_{j}}\\
\frac{\partial f_{2}}{\partial p_{j}}\\
\vdots\\
\frac{\partial f_{k}}{\partial p_{j}}
\end{array}\right)
```

are the parameter derivatives, and

```math
S_{j}=\left(\begin{array}{c}
\frac{\partial y_{1}}{\partial p_{j}}\\
\frac{\partial y_{2}}{\partial p_{j}}\\
\vdots\\
\frac{\partial y_{k}}{\partial p_{j}}
\end{array}\right)
```

is the vector of sensitivities. Since this ODE is dependent on the values of the
independent variables themselves, this ODE is computed simultaneously with the
actual ODE system.

Note that the Jacobian-vector product:

```math
\frac{\partial f}{\partial u}\frac{\partial u}{\partial p_{j}}
```

can be computed without forming the Jacobian. With finite differences, this through using the following
formula for the directional derivative:

```math
Jv \approx \frac{f(x+v \epsilon) - f(x)}{\epsilon}
```

or by using a dual number with a single partial dimension, ``d = x + v \epsilon`` we get that

```math
f(d) = f(x) + Jv \epsilon
```

as a fast way to calcuate ``Jv``. Thus, except when a sufficiently good function for `J` is given
by the user, the Jacobian is never formed. For more details, consult the 
[MIT 18.337 lecture notes on forward mode AD](https://mitmath.github.io/18337/lecture9/autodiff_dimensions).

#### Syntax

`ODELocalSensitivityProblem` is similar to an `ODEProblem`:

```julia
function ODELocalSensitivityProblem(f::DiffEqBase.AbstractODEFunction,u0,
                                    tspan,p=nothing,
                                    SensitivityAlg();
                                    kwargs...)
```

The `SensitivityAlg` is used to mirror the definition of the derivative options seen generally
throughout OrdinaryDiffEq.jl. The keyword options on the `SensitivityAlg` are as follows:

* `autojacvec`: Calculate Jacobian-vector (local sensitivity analysis) product 
  via automatic differentiation with special seeding to avoid building the Jacobian. 
  Default is `true`. If `autodiff=false`, it will use the Jacobian-free forward
  differencing approximation. If `false`, the Jacobian will prefer to use any 
  `f.jac` function provided by the user. If none is provided by the user, then it 
  will fall back to automatic or finite differentiation, though this choice is 
  not recommended.
* `autodiff`: Use automatic differentiation in the internal sensitivity algorithm
  computations. Default is `true`.
* `chunk_size`: Chunk size for forward mode differentiation. Default is `0` for
  automatic chunk size choice. Only used when `autojacvec=false`.
* `diff_type`: Choice of differencing used to build the Jacobian when `autojacvec=false`
  and `autodiff=false`. Defaults to `Val{:central}` for central differencing with
  DiffEqDiffTools.jl.

#### Example solving an ODELocalSensitivityProblem

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
where `dp[i]` is the gradient of component `i` by the parameters. The first gives
the full timeseries of values. The second returns the `i`th values, while the third
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
autodifferentiation and numerical differentiation through the ODE solver:

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

For completeness, we detail the internal representation. Therefore, the solution
to the ODE are the first `n` components of the solution. This means we can grab
the matrix of solution values like:

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

## Adjoint Sensitivity Analysis

Adjoint sensitivity analysis is used to find the gradient of the solution
with respect to some functional of the solution. In many cases this is used
in an optimization problem to return the gradient with respect to some cost
function. It is equivalent to "backpropogation" or reverse-mode automatic
differentiation of a differential equation.

### Adjoint Sensitivity Analysis via adjoint_sensitivities

This adjoint requires the definition of some scalar functional ``g(u,p,t)``
where ``u`` is the (numerical) solution to the differential equation.
Adjoint sensitivity analysis finds the gradient of

```math
G(u,p)=G(u(p))=\int_{t_{0}}^{T}g(u(t,p))dt
```

some integral of the solution. It does so by solving the adjoint problem

```math
\frac{d\lambda^{\star}}{dt}=g_{u}(t)-\lambda^{\star}(t)f_{u}(t),\thinspace\thinspace\thinspace\lambda^{\star}(T)=0
```

where ``f_u`` is the Jacobian of the system with respect to the state `u` while
``f_p`` is the Jacobian with respect to the parameters. The adjoint problem's
solution gives the sensitivities through the integral:

```math
\frac{dG}{dp}=\int_{t_{0}}^{T}\lambda^{\star}(t)f_{p}(t)+g_{p}(t)dt+\lambda^{\star}(t_{0})u_{p}(t_{0})
```

Notice that since the adjoints require the Jacobian of the system at the state,
it requires the ability to evaluate the state at any point in time. Thus it
requires the continuous forward solution in order to solve the adjoint solution,
and the adjoint solution is required to be continuous in order to calculate the
resulting integral.

There is one extra detail to consider. In many cases we would like to calculate
the adjoint sensitivity of some discontinuous functional of the solution. One
canonical function is the L2 loss against some data points, that is:

```math
L(u,p,t)=\sum_{i=1}^{n}\Vert\tilde{u}(t_{i})-u(t_{i},p)\Vert^{2}
```

In this case, we can reinterpret our summation as the distribution integral:

```math
G(u,p)=\int_{0}^{T}\sum_{i=1}^{n}\Vert\tilde{u}(t_{i})-u(t_{i},p)\Vert^{2}\delta(t_{i}-t)dt
```

where ``δ`` is the Dirac distribution. In this case, the integral is continuous
except at finitely many points. Thus it can be calculated between each ``t_i``.
At a given ``t_i``, given that the ``t_i`` are unique, we have that

```math
g_{y}(t_{i})=2\left(\tilde{u}(t_{i})-u(t_{i},p)\right)
```

Thus the adjoint solution is given by integrating between the integrals and
applying the jump function ``g_y`` at every data point.

We note that

```math
\lambda^{\star}(t)f_{u}(t)
```

is a vector-transpose Jacobian product, also known as a `vjp`, which can be efficiently computed 
using the pullback of backpropogation on the user function `f` with a forward pass at `u` with a
pullback vector ``\lambda^{\star}``. For more information, consult the
[MIT 18.337 lecture notes on reverse mode AD](https://mitmath.github.io/18337/lecture10/estimation_identification)

#### Syntax

There are two forms. For discrete adjoints, the form is:

```julia
s = adjoint_sensitivities(sol,alg,dg,ts;kwargs...)
```

where `alg` is the ODE algorithm to solve the adjoint problem, `dg` is the jump
function, and `ts` is the time points for data. `dg` is given by:

```julia
dg(out,u,p,t,i)
```

which is the in-place gradient of the cost functional `g` at time point `ts[i]`
with `u=u(t)`.

For continuous functionals, the form is:

```julia
s = adjoint_sensitivities(sol,alg,g,nothing,dg;kwargs...)
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
s = adjoint_sensitivities(sol,alg,g,nothing;kwargs...)
```

then it will be computed automatically using ForwardDiff or finite
differencing, depending on the `autodiff` setting in the `SensitivityAlg`.
Note that the keyword arguments are passed to the internal ODE solver for 
solving the adjoint problem. Two special keyword arguments are `iabstol` and 
`ireltol` which are the tolerances for the internal quadrature via QuadGK for 
the resulting functional.

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

#### Example discrete adjoints on a cost function

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
\begin{align}
dg_{1}&=1-u_{1} \\
dg_{2}&=1-u_{2}
\end{align}
```

and thus:

```julia
dg(out,u,i) = (out.=1.0.-u)
```

If we had data, we'd just replace `1.0` with `data[i]`. To get the adjoint
sensitivities, call:

```julia
res = adjoint_sensitivities(sol,Vern9(),dg,t,abstol=1e-14,
                            reltol=1e-14,iabstol=1e-14,ireltol=1e-12)
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

#### Example controlling adjoint method choices and checkpointing

In the previous examples, all calculations were done using the interpolating
method. This maximizes speed but at a cost of requiring a dense `sol`. If it
is not possible to hold a dense forward solution in memory, then one can use
checkpointing. This is enabled by default if `sol` is not dense, so for
example

```julia
sol = solve(prob,Vern9(),saveat=[0.0,0.2,0.5,0.7])
```

Creates a non-dense solution with checkpoints at `[0.0,0.2,0.5,0.7]`. Now we
can do

```julia
res = adjoint_sensitivities(sol,Vern9(),dg,t)
```

When grabbing a Jacobian value during the backwards solution, it will no longer
interpolate to get the value. Instead, it will start a forward solution at the
nearest checkpoint and solve until the necessary time.

To eliminate the extra forward solutions, one can instead pass the `SensitivityAlg`
with the `backsolve=true` option:

```julia
sol = solve(prob,Vern9(),save_everystep=false,save_start=false)
res = adjoint_sensitivities(sol,Vern9(),dg,t,sensealg=SensitivityAlg(backsolve=true))
```

When this is done, the values for the Jacobian will be computing the original ODE in
reverse. Note that this only requires the final value of the solution.

#### Applicability of Backsolve and Caution

When `backsolve` is applicable it is the fastest method and requires the least memory.
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
enabling this method.

#### Example continuous adjoints on an energy functional

In this case we'd like to calculate the adjoint sensitivity of the scalar energy
functional

```math
G(u,p)=\int_{0}^{T}\frac{\sum_{i=1}^{n}u_{i}^{2}(t)}{2}dt
```

which is

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
