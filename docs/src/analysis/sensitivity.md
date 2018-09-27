# Sensitivity Analysis

Sensitivity analysis for ODE models is provided by the DiffEq suite. The model
sensitivities are defined as the derivatives of the solution with respect to the
parameters. Sensitivity analysis serves two major purposes. On one hand, the
sensitivities are diagnostics of the model which are useful for understand how
it will change in accordance to changes in the parameters. But another use is
simply because in many cases these derivatives are useful. Sensitivity analysis
provides a cheap way to calculate the gradient of the solution which can be
used in parameter estimation and other optimization tasks.

There are three types of sensitivity analysis. Local sensitivity analysis directly
gives the gradient of the solution with respect to each parameter along the time
series. The computational cost scales like `N*M`, where `N` is the number of states
and `M` is the number of parameters. While this gives all of the information,
it can be expensive for large models. Instead, adjoint sensitivity analysis solves
directly for the gradient of some functional of the solution, such as a cost
function or energy functional, in a much cheaper manner. Global Sensitivty Analysis methods
are meant to be used for exploring the sensitivity over a larger domain without calculating
derivatives.

#### Note

Currently there are more performance optimizations needed to be done on the
adjoint sensitivity method.

## Local Sensitivity Analysis

The local sensitivity of the solution to a parameter is defined by how much the
solution would change by changes in the parameter, i.e. the sensitivity of the
ith independent variable to the jth parameter is `` \frac{\partial u}{\partial p_{j}}``.

The local sensitivity is computed using the sensitivity ODE:

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

### Example solving an ODELocalSensitivityProblem

To define a sensitivity problem, simply use the `ODELocalSensitivityProblem` type
instead of an ODE type. Note that this requires a [ParameterizedFunction](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl) with a
Jacobian. For example, we generate an ODE with the sensitivity equations attached
for the Lotka-Volterra equations by:

```julia
f = @ode_def_nohes LotkaVolterraSensitivity begin
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
function.

The adjoint requires the definition of some scalar functional ``g(u,p,t)``
where ``u`` is the (numerical) solution to the differential equation.
Adjoint sensitivity analysis finds the gradient of

```math
G(u,p)=G(u(p))=\int_{t_{0}}^{T}g(u(t,p))dt
```

some integral of the solution. It does so by solving the adjoint problem

```math
\frac{d\lambda^{\star}}{dt}=-g_{u}(t)-\lambda^{\star}(t)f_{u}(t),\thinspace\thinspace\thinspace\lambda^{\star}(T)=0
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

### Syntax

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

### Example discrete adjoints on a cost function

In this example we will show solving for the adjoint sensitivities of a discrete
cost functional. First let's solve the ODE and get a high quality continuous
solution:

```julia
f = @ode_def_nohes LotkaVolterra begin
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

## Global Sensitivity Analysis

Global Sensitivity Analysis methods are used to quantify the uncertainity in
output of a model w.r.t. the parameters, their individual contributions or the
contribution of their interactions. The type of GSA method to use depends on
the interest of the user, below we describe the methods available in the suite
at the moment (some more are already in development) and explain what is
the output of each of the methods and what it represents.

### Morris Method

The Morris method also known as Morris’s OAT method where OAT stands for
One At a Time can be described in the following steps:

We calculate local sensitivity measures known as “elementary effects”,
which are calculated by measuring the perturbation in the output of the
model on changing one parameter.

``EE_i = \frac{f(x_1,x_2,..x_i+ \Delta,..x_k) - y}{\Delta}``

These are evaluated at various points in the input chosen such that a wide
“spread” of the parameter space is explored and considered in the analysis,
to provide an approximate global importance measure. The mean and variance of
these elementary effects is computed. A high value of the mean implies that
a parameter is important, a high variance implies that its effects are
non-linear or the result of interactions with other inputs. This method
does not evaluate separately the contribution from the
interaction and the contribution of the parameters individually and gives the
effects for each parameter which takes into cpnsideration all the interactions and its
individual contribution.

`morris_effects = morris_sensitivity(f,param_range,param_steps;relative_scale=false,kwargs...)`

`morris_effects = morris_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,param_steps;kwargs...)`

Here, `f` is just the model (as a julia function or a `DEProblem`) you want to
run the analysis on, `param_range` requires an array of 2-tuples with the lower bound
and the upper bound, `param_steps` decides the value of \Delta in the equation
above and `relative_scale`, the above equation takes the assumption that
the parameters lie in the range `[0,1]` but as this is not always the case
scaling is used to get more informative, scaled effects.

### Sobol Method

Sobol is a variance-based method and it decomposes the variance of the output of
the model or system into fractions which can be attributed to inputs or sets
of inputs. This helps to get not just the individual parameter's sensitivities
but also gives a way to quantify the affect and sensitivity from
the interaction between the parameters.

```math
 Y = f_0+ \sum_{i=1}^d f_i(X_i)+ \sum_{i < j}^d f_{ij}(X_i,X_j) ... + f_{1,2...d}(X_1,X_2,..X_d)
```

```math
 Var(Y) = \sum_{i=1}^d V_i + \sum_{i < j}^d V_{ij} + ... + V_{1,2...,d}
```

The Sobol Indices are "order"ed, the first order indices given by ``S_i = \frac{V_i}{Var(Y)}``
the contribution to the output variance of the main effect of `` X_i ``, therefore it
measures the effect of varying `` X_i `` alone, but averaged over variations
in other input parameters. It is standardised by the total variance to provide a fractional contribution.
Higher-order interaction indices `` S_{i,j}, S_{i,j,k} `` and so on can be formed
by dividing other terms in the variance decomposition by `` Var(Y) ``.

`sobol_second_order = sobol_sensitivity(f,param_range,N,order=2)`

`sobol_second_order = sobol_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,N,order=2)`

Here `f` and `param_range` are the same as Morris's, providing a uniform interface.

### Regression Method

If a sample of inputs and outputs `` (X^n, Y^n) = 􏰀(X^{i}_1, . . . , X^{i}_d, Y_i)_{i=1..n} ``􏰁
is available, it is possible to fit a linear model explaining the behaviour of Y given the
values of X, provided that the sample size n is sufficiently large (at least n > d).

The measures provided for this analysis by us in DiffEqSensitivity.jl are

  a) Pearson Correlation Coefficient:

```math
r = \frac{\sum_{i=1}^{n} (x_i - \overline{x})(y_i - \overline{y})}{\sqrt{\sum_{i=1}^{n} (x_i - \overline{x})^2(y_i - \overline{y})^2}}
```

  b) Standard Regression Coefficient (SRC):

```math
SRC_j = \beta_{j} \sqrt{\frac{Var(X_j)}{Var(Y)}}
```

  where ``\beta_j`` is the linear regression coefficient associated to $X_j$.

  c) Partial Correlation Coefficient (PCC):

```math
PCC_j = \rho(X_j - \hat{X_{-j}},Y_j - \hat{Y_{-j}})
```

  where ``\hat{X_{-j}}``􏰈 is the prediction of the linear model, expressing ``X_{j}``
  with respect to the other inputs and ``\hat{Y􏰈_{-j}}`` is the prediction of the
  linear model where ``X_j`` is absent. PCC measures the sensitivity of ``Y`` to
  ``X_j`` when the effects of the other inputs have been canceled.

`regre_sensitivity = regression_sensitivity(f,param_range,param_fixed,n;coeffs=:rank)`

`regre_sensitivity = regression_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,param_fixed,n;coeffs=:rank)`

Again, `f` and `param_range` are the same as above. An array of the true parameter values
that lie within the `param_range` bounds are passed through the `param_fixed` argument.
`n` determines the number of simulations of the model run to generate the data points
of the solution and parameter values and the `coeffs` kwarg lets you decide the
coefficients you want.

### GSA example

Let's create the ODE problem to run our GSA on.

```julia
f = @ode_def_nohes LotkaVolterraTest begin
    dx = a*x - b*x*y
    dy = -3*y + x*y
end a b
u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5,1.0]
prob = ODEProblem(f,u0,tspan,p)
t = collect(range(0, stop=10, length=200))
```
For Morris Method

```julia
m = DiffEqSensitivity.morris_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],[10,10],len_trajectory=1500,total_num_trajectory=1000,num_trajectory=150)
```
Let's get the means and variances from the `MorrisSensitivity` struct.

```julia
m.means

Out[9]: 2-element Array{Array{Float64,2},1}:
 [0.0 0.0513678 … 7.91336 7.93783; 0.0 0.00115769 … 3.66156 3.67284]
 [0.0 0.0488899 … 2.50728 2.359; 0.0 0.00112006 … 2.23431 2.44946]

m.variances

Out[10]: 2-element Array{Array{Float64,2},1}:
 [0.0 1.94672e-5 … 26.4223 24.8513; 0.0 4.81347e-9 … 37.4061 30.3068]
 [0.0 1.77615e-5 … 17.9555 14.9231; 0.0 4.47931e-9 … 48.074 51.9312]
```
This gives the means of the effects and it's variances over the entire timespan and thus we get 200-length
arrays for each paramter and dependent variable pair.

We can plot the trajectory of the sensitivity with the standard deviation bars.
```julia
# For the first parameter (a)
stdv1 = sqrt.(m.variances[1])
p = plot(m.means[1]', yerror=stdv1)
```
![morrisparameter1](../assets/morris1.png)

```julia
# For the second parameter (b)
stdv2 = sqrt.(m.variances[2])
p = plot(m.means[2]', yerror=stdv2)
```
![morrisparameter2](../assets/morris2.png)

For Sobol Method

```julia

s0 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,0)
Out[8]: 2-element Array{Array{Float64,2},1}:
 [NaN 0.507831 … 1.00731 1.00436; NaN 1.92336 … 0.732384 0.730945]
 [NaN 0.47214 … 0.676224 0.681525; NaN -1.68656 … 0.879557 0.877603]

s1 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,1)
Out[9]: 2-element Array{Array{Float64,2},1}:
 [NaN 0.39537 … 0.341697 0.343645; NaN -2.06101 … 0.10922 0.106976]
 [NaN 0.652815 … 0.00910675 0.00815206; NaN 5.24832 … 0.296978 0.296639]

s2 = sobol_sensitivity(prob,Tsit5(),t,[[1,5],[0.5,5]],N,2)
Out[10]: 1-element Array{Array{Float64,2},1}:
 [NaN -0.0596478 … 0.652303 0.657847; NaN -1.84504 … 0.645139 0.620036]
```
We can decide which order of Sobol Indices we are interested in my passing an argument for it,
by default it gives the second order indices. Again the result is obtained over the entire `timespan`

We plot the first order and total order Sobol Indices for some timepoints for each of the parameters (`a` and `b`).

```julia

p1 = bar(["a","b"],[s0[1][end-2],s0[2][end-2]],color=[:red,:blue],title="Total Order Indices at t=9.949748743718592",legend=false)
p2 = bar(["a","b"],[s1[1][end-2],s1[2][end-2]],color=[:red,:blue],title="First Order Indices at t=9.949748743718592",legend=false)
p3 = bar(["a","b"],[s0[1][3],s0[2][3]],color=[:red,:blue],title="Total Order Indices at t=0.05025125628140704",legend=false)
p4 = bar(["a","b"],[s1[1][3],s1[2][3]],color=[:red,:blue],title="First Order Indices at t=0.05025125628140704",legend=false)
plo = plot(p1,p2,p3,p4,layout=(4,1),size=(600,500))

```
![sobolplot](../assets/sobolbars.png)

Here we plot the Sobol indices of first order and the total Sobol indices for the parameters `a` and `b`. The plots are obtained by getting the Sobol Indices at the `t = 9.949748743718592` and the `t = 0.05025125628140704` time point of the first dependent variable `x(t)` from the 200-length sensitivities over the entire time span. The length of the bar represents the quantification of the sensitivity of the output to that parameter and hence for the 199th time point you can say that `x(t)` is more sensitive to `b`, also you can observe how the relative difference between `a` and `b` is larger in the first order than the total order indices, this tells us that most of the contribution of `a` to `x(t)` arises from interactions and it's individual non-interaction contribution is significantly lesser than `b` and vice-versa for `b` as it's first order plot indicates quite high value.

