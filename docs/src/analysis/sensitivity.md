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

## Efficiency of the Different Methods

For an analysis of which methods will be most efficient for computing the
solution derivatives for a given problem, consult our analysis
[in this arxiv paper](https://arxiv.org/abs/1812.01892). A general rule of thumb
is:

- Discrete Forward Sensitivity Analysis via ForwardDiff.jl is the fastest for
  ODEs with small numbers of parameters (<100)
- Adjoint senstivity analysis is the fastest when the number of parameters is
  sufficiently large.
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
f = @ode_def begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c

p = [1.5,1.0,3.0,1.0]
u0 = [1.0;1.0]
prob = ODEProblem(f,u0,(0.0,10.0),p)
```

Let's say we wanted to get the derivative of the final value w.r.t. each of
the parameters. We can define the following function:

```julia
function test_f(p)
  _prob = remake(prob;u0=convert.(eltype(p),prob.u0),p=p)
  solve(prob,Vern9(),save_everystep=false)[end]
end
```

What this function does is use the `remake` function
[from the Problem Interface page](http://docs.juliadiffeq.org/latest/basics/problem.html#Modification-of-problem-types-1)
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
p4dual = Dual{MyTag}(3.0, (0.0, 0.0, 0.0, 0.0))
pdual = [p1dual, p2dual, p3dual, p4dual]
```

or equivalently using the `seed_duals` convenience function:

```julia
pdual = seed_duals(p,MyTag)
```

Next we need to make our initial condition Dual numbers so that these propogate
through the solution. We can do this manually like:

```julia
u0dual = [Dual{Float64}(1.0, (0.0, 0.0)),Dual{Float64}(1.0, (0.0, 0.0))]
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

#### Example solving an ODELocalSensitivityProblem

To define a sensitivity problem, simply use the `ODELocalSensitivityProblem` type
instead of an ODE type. For example, we generate an ODE with the sensitivity
equations attached for the Lotka-Volterra equations by:

```julia
f = @ode_def begin
  dx = a*x - b*x*y
  dy = -c*y + x*y
end a b c

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

Currently, the gradient is required. Note that the keyword arguments are passed
to the internal ODE solver for solving the adjoint problem. Two special keyword
arguments are `iabstol` and `ireltol` which are the tolerances for the internal
quadrature via QuadGK for the resulting functional.

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

* `quad`: Use Gauss-Kronrod quadrature to integrate the adjoint sensitivity
  integral. Disabling it can decrease memory usage but increase computation
  time. Default is `true`.
* `backsolve`: Solve the differential equation backward to get the past values.
  Note that for chaotic or non-reversible systems, enabling this option can
  lead to wildly incorrect results. Enabling it can decrease memory usage but
  increase computation time. When it is set to `true`, `quad` will be
  automatically set to `false`. Default is `false`.
* `autodiff`: Use automatic differentiation. Default is `true`.
* `chunk_size`: Chunk size for forward mode differentiation. Default is `0`.
* `autojacvec`: Calculate Jacobian-vector (local sensitivity analysis) or
  vector-Jacobian (adjoint sensitivity analysis) product via automatic
  differentiation with special seeding. Default is `true`.

#### Example discrete adjoints on a cost function

In this example we will show solving for the adjoint sensitivities of a discrete
cost functional. First let's solve the ODE and get a high quality continuous
solution:

```julia
f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + x*y
end a b c

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
  tmp_prob = problem_new_parameters(prob,p)
  sol = solve(tmp_prob,Vern9(),abstol=1e-14,reltol=1e-14)
  res,err = quadgk((t)-> (sum(sol(t)).^2)./2,0.0,10.0,abstol=1e-14,reltol=1e-10)
  res
end
res2 = ForwardDiff.gradient(G,[1.5,1.0,3.0])
res3 = Calculus.gradient(G,[1.5,1.0,3.0])
```
