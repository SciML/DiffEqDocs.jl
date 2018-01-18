# Sensitivity Analysis

Sensitivity analysis for ODE models is provided by the DiffEq suite. The model
sensitivities are defined as the derivatives of the solution with respect to the
parameters. Sensitivity analysis serves two major purposes. On one hand, the
sensitivities are diagnostics of the model which are useful for understand how
it will change in accordance to changes in the parameters. But another use is
simply because in many cases these derivatives are useful. Sensitivity analysis
provides a cheap way to calculate the gradient of the solution which can be
used in parameter estimation and other optimization tasks.

There are two types of sensitivity analysis. Local sensitivity analysis directly
gives the gradient of the solution with respect to each parameter along the time
series. The computational cost scales like `N*M`, where `N` is the number of states
and `M` is the number of parameters. While this gives all of the information,
it can be expensive for large models. Instead, adjoint sensitivity analysis solves
directly for the gradient of some functional of the solution, such as a cost
function or energy functional, in a much cheaper manner.

#### Note

Currently there are more performance optimizations needed to be done on the
adjoint sensitivity method.

## Local Sensitivity Analysis

The local sensitivity of the solution to a parameter is defined by how much the
solution would change by changes in the parameter, i.e. the sensitivity of the
ith independent variable to the jth parameter is `` \frac{\partial y}{\partial p_{j}}``.

The local sensitivity is computed using the sensitivity ODE:

```math
\frac{d}{dt}\frac{\partial u}{\partial p_{j}}=\frac{\partial f}{\partial y}\frac{\partial y}{\partial p_{j}}+\frac{\partial f}{\partial p_{j}}=J\cdot S_{j}+F_{j}
```

where

```math
J=\left(\begin{array}{cccc}
\frac{\partial f_{1}}{\partial y_{1}} & \frac{\partial f_{1}}{\partial y_{2}} & \cdots & \frac{\partial f_{1}}{\partial y_{k}}\\
\frac{\partial f_{2}}{\partial y_{1}} & \frac{\partial f_{2}}{\partial y_{2}} & \cdots & \frac{\partial f_{2}}{\partial y_{k}}\\
\cdots & \cdots & \cdots & \cdots\\
\frac{\partial f_{k}}{\partial y_{1}} & \frac{\partial f_{k}}{\partial y_{2}} & \cdots & \frac{\partial f_{k}}{\partial y_{k}}
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

### Example solving an ODELocalSensitivityProblem

To define a sensitivity problem, simply use the `ODELocalSensitivityProblem` type
instead of an ODE type. Note that this requires a [ParameterizedFunction](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl) with a
Jacobian. For example, we generate an ODE with the sensitivity equations attached
for the Lotka-Volterra equations by:

```julia
f = @ode_def_nohes LotkaVolterraSensitivity begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=>3 d=1

prob = ODELocalSensitivityProblem(f,[1.0;1.0],(0.0,10.0))
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
  pf = ParameterizedFunction(f,p)
  prob = ODEProblem(pf,eltype(p).([1.0,1.0]),eltype(p).((0.0,10.0)))
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
function.

The adjoint requires the definition of some scalar functional ``g(t,u,p)``
where ``u`` is the (numerical) solution to the differential equation.
Adjoint sensitivity analysis finds the gradient of

```math
G(u,p)=G(u(p))=\int_{0}^{T}g(u(t,p))dt
```

some integral of the solution. It does so by solving the adjoint problem

```math
\frac{d\lambda^{\star}}{dt}=g_{u}(t)-\lambda^{\star}(t)f_{u}(t),\thinspace\thinspace\thinspace\lambda^{\star}=0
```

where ``f_u`` is the Jacobian of the system with respect to the state `u` while
``f_p`` is the Jacobian with respect to the parameters. The adjoint problem's
solution gives the sensitivities through the integral:

```math
\frac{dG}{dp}=-\int_{0}^{T}\lambda^{\star}(t)f_{p}(t)dt-\lambda^{T}(0)u_{p}(t)
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
L(t,u,p)=\sum_{i=1}^{n}\Vert\tilde{u}(t_{i})-u(t_{i},p)\Vert^{2}
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

### Syntax

There are two forms. For discrete adjoints, the form is:

```julia
s = adjoint_sensitivities(sol,alg,dg,t;kwargs...)
```

where `alg` is the ODE algorithm to solve the adjoint problem, `dg` is the jump
function, and `t` is the time points for data. `dg` is given by:

```julia
dg(out,u,i)
```

which is the in-place gradient of the cost functional `g` at time point `t[i]`
with `u=u(t[i])`.

For continuous functionals, the form is:

```julia
s = adjoint_sensitivities(sol,alg,g,nothing,dg;kwargs...)
```

for the cost functional

```julia
g(t,u,p)
```

with in-place gradient

```julia
dg(out,t,u,p)
```

Currently, the gradient is required. Note that the keyword arguments are passed
to the internal ODE solver for solving the adjoint problem. Two special keyword
arguments are `iabstol` and `ireltol` which are the tolerances for the internal
quadrature via QuadGK for the resulting functional.

### Example discrete adjoints on a cost function

In this example we will show solving for the adjoint sensitivities of a discrete
cost functional. First let's solve the ODE and get a high quality continuous
solution:

```julia
f = @ode_def_nohes LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1.0 c=>3.0 d=1

prob = ODEProblem(f,[1.0;1.0],(0.0,10.0))
sol = solve(prob,Vern9(),abstol=1e-10,reltol=1e-10)
```

Now let's calculate the sensitivity of the L2 error against 1 at evenly spaced
points in time, that is:

```math
L(t,u,p)=\sum_{i=1}^{n}\frac{\Vert1-u(t_{i},p)\Vert^{2}}{2}
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
  tmp_prob = problem_new_parameters(prob,p)
  sol = solve(tmp_prob,Vern9(),abstol=1e-14,reltol=1e-14,saveat=t)
  A = convert(Array,sol)
  sum(((1-A).^2)./2)
end
G([1.5,1.0,3.0])
res2 = ForwardDiff.gradient(G,[1.5,1.0,3.0])
res3 = Calculus.gradient(G,[1.5,1.0,3.0])
```

and see this gives the same values.

### Example continuous adjoints on an energy functional

In this case we'd like to calculate the adjoint sensitivity of the scalar energy
functional

```math
G(u,p)=\int_{0}^{T}\frac{\sum_{i=1}^{n}u_{i}^{2}(t)}{2}dt
```

which is

```julia
g(t,u,p) = (sum(u).^2) ./ 2
```

Notice that the gradient of this function with respect to the state `u` is:

```julia
function dg(out,t,u,p)
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
