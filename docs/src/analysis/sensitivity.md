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

$$ EE_i = \frac{f(x_1,x_2,..x_i+ \Delta,..x_k) - y}{\Delta} $$

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
run the analysis on, `param_range` requires an array of your range for the 
parameters as an array of the lower bound and the upper bound, `param_steps` decides the value of \Delta in the equation 
above and `relative_scale`, the above equation takes the assumption that 
the parameters lie in the range `[0,1]` but as this is not always the case 
scaling is used to get more informative, scaled effects.

### Sobol Method

Sobol is a variance-based method and it decomposes the variance of the output of 
the model or system into fractions which can be attributed to inputs or sets 
of inputs. This helps to get not just the individual parameter's sensitivities 
but also gives a way to quantify the affect and sensitivity from 
the interaction between the parameters. 

$$ Y = f_0+ \sum_{i=1}^d f_i(X_i)+ \sum_{i < j}^d f_{ij}(X_i,X_j) ... + f_{1,2...d}(X_1,X_2,..X_d) $$
                                   
$$ Var(Y) = \sum_{i=1}^d V_i + \sum_{i < j}^d V_{ij} + ... + V_{1,2...,d} $$

The Sobol Indices are "order"ed, the first order indices given by $$S_i = \frac{V_i}{Var(Y)}$$ 
the contribution to the output variance of the main effect of $$ X_i $$, therefore it 
measures the effect of varying $$ X_i $$ alone, but averaged over variations 
in other input parameters. It is standardised by the total variance to provide a fractional contribution. 
Higher-order interaction indices $$ S_{i,j}, S_{i,j,k} $$ and so on can be formed 
by dividing other terms in the variance decomposition by $$ Var(Y) $$.

`sobol_second_order = sobol_sensitivity(f,param_range,N,order=2)`

`sobol_second_order = sobol_sensitivity(prob::DiffEqBase.DEProblem,alg,t,param_range,N,order=2)`

Here `f` and `param_range` are the same as Morris's, providing a uniform interface.

### Regression Method

If a sample of inputs and outputs $$ (X^n, Y^n) = 􏰀(X^{i}_1, . . . , X^{i}_d, Y_i)_{i=1..n} $$􏰁 
is available, it is possible to fit a linear model explaining the behaviour of Y given the 
values of X, provided that the sample size n is sufficiently large (at least n > d).

The measures provided for this analysis by us in DiffEqSensitivity.jl are

  a) Pearson Correlation Coefficient:

$$ r = \frac{\sum_{i=1}^{n} (x_i - \overline{x})(y_i - \overline{y})}
{\sqrt{\sum_{i=1}^{n} (x_i - \overline{x})^2(y_i - \overline{y})^2}} $$

  b) Standard Regression Coefficient (SRC):

$$ SRC_j = \beta_{j} \sqrt{\frac{Var(X_j)}{Var(Y)}} $$

  where $$ \beta_j $$ is the linear regression coefficient associated to $X_j$.

  c) Partial Correlation Coefficient (PCC):

$$ PCC_j = \rho(X_j - \hat{X_{-j}},Y_j - \hat{Y_{-j}}) $$

  where $$ \hat{X_{-j}} $$􏰈 is the prediction of the linear model, expressing $$ X_{j} $$ 
  with respect to the other inputs and $$ \hat{Y􏰈_{-j}} $$ is the prediction of the 
  linear model where $$ X_j $$ is absent. PCC measures the sensitivity of $$ Y $$ to 
  $$ X_j $$ when the effects of the other inputs have been canceled.

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


Out[8]:MorrisSensitivity(Array{Float64,2}[[0.0 0.0513678 … 7.91336 7.93783; 0.0 0.00115769 … 3.66156 3.67284], [0.0 0.0488899 … 2.50728 2.359; 0.0 0.00112006 … 2.23431 2.44946]], Array{Float64,2}[[0.0 1.94672e-5 … 26.4223 24.8513; 0.0 4.81347e-9 … 37.4061 30.3068], [0.0 1.77615e-5 … 17.9555 14.9231; 0.0 4.47931e-9 … 48.074 51.9312]], Array{Array{Float64,2},1}[[[0.0 0.0502074 … 7.98867 9.97645; 0.0 0.00113922 … 2.91223 2.50633], [0.0 0.0580075 … 11.2991 7.40316; 0.0 0.0012612 … 5.7219 7.9419], [0.0 0.0552782 … 11.1432 12.8587; 0.0 0.00121893 … 5.19252 4.16375], [0.0 0.0540595 … 11.0023 11.2203; 0.0 0.0012001 … 1.14073 1.73876], [0.0 0.0541561 … 1.73364 2.29789; 0.0 0.00120179 … 2.46675 2.15942], [0.0 0.0541561 … 1.73364 2.29789; 0.0 0.00120179 … 2.46675 2.15942], [0.0 0.0554766 … 0.664665 1.09787; 0.0 0.00122238 … 3.62296 3.2128], [0.0 0.0554766 … 0.664665 1.09787; 0.0 0.00122238 … 3.62296 3.2128], [0.0 0.0554766 … 0.664665 1.09787; 0.0 0.00122238 … 3.62296 3.2128], [0.0 0.0541561 … 1.73364 2.29789; 0.0 0.00120179 … 2.46675 2.15942]  …  [0.0 0.0595316 … 0.322627 0.135043; 0.0 0.00128482 … 9.60211 8.46034], [0.0 0.0581124 … 1.99996 2.50206; 0.0 0.001263 … 2.12335 1.79987], [0.0 0.0608753 … 0.532116 0.708165; 0.0 0.00130521 … 3.34051 2.87788], [0.0 0.0608753 … 0.532116 0.708165; 0.0 0.00130521 … 3.34051 2.87788], [0.0 0.0609864 … 20.0782 13.8505; 0.0 0.00130709 … 20.0936 37.0211], [0.0 0.0609864 … 20.0782 13.8505; 0.0 0.00130709 … 20.0936 37.0211], [0.0 0.0623634 … 30.6292 28.7567; 0.0 0.00132789 … 65.3938 39.8315], [0.0 0.0623634 … 30.6292 28.7567; 0.0 0.00132789 … 65.3938 39.8315], [0.0 0.0608753 … 0.532116 0.708165; 0.0 0.00130521 … 3.34051 2.87788], [0.0 0.0608753 … 0.532116 0.708165; 0.0 0.00130521 … 3.34051 2.87788]], [[0.0 0.0488908 … 0.559981 0.41625; 0.0 0.00111992 … 0.309004 0.245778], [0.0 0.0478053 … 1.13904 0.800567; 0.0 0.00110272 … 0.371875 0.471092], [0.0 0.0489784 … 0.777046 0.447304; 0.0 0.0011215 … 0.409923 0.519859], [0.0 0.0501811 … 0.509473 0.491546; 0.0 0.00114065 … 0.26971 0.357191], [0.0 0.0514142 … 0.0317235 0.388981; 0.0 0.00116017 … 0.0443424 0.0666987], [0.0 0.0526785 … 1.34268 0.635881; 0.0 0.00118008 … 0.381449 0.521107], [0.0 0.0539748 … 4.35505 4.39203; 0.0 0.0012004 … 0.409809 0.657726], [0.0 0.0552028 … 11.2465 11.0146; 0.0 0.00121937 … 2.61037 0.503488], [0.0 0.0538765 … 3.00441 1.96804; 0.0 0.00119868 … 2.28968 2.68301], [0.0 0.0515072 … 0.539924 0.601361; 0.0 0.00116181 … 0.191549 0.165097]  …  [0.0 0.0566669 … 2.35688 2.84764; 0.0 0.00124227 … 4.30777 3.66766], [0.0 0.0554054 … 0.292358 0.503547; 0.0 0.00122288 … 10.9555 9.58808], [0.0 0.0566669 … 2.35688 2.84764; 0.0 0.00124227 … 4.30777 3.66766], [0.0 0.0580645 … 29.528 26.3741; 0.0 0.00126386 … 13.7705 29.949], [0.0 0.0580645 … 29.528 26.3741; 0.0 0.00126386 … 13.7705 29.949], [0.0 0.0580645 … 29.528 26.3741; 0.0 0.00126386 … 13.7705 29.949], [0.0 0.0580645 … 29.528 26.3741; 0.0 0.00126386 … 13.7705 29.949], [0.0 0.0580645 … 29.528 26.3741; 0.0 0.00126386 … 13.7705 29.949], [0.0 0.0580645 … 29.528 26.3741; 0.0 0.00126386 … 13.7705 29.949], [0.0 0.0566669 … 2.35688 2.84764; 0.0 0.00124227 … 4.30777 3.66766]]])
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



