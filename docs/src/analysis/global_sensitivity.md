# Global Sensitivity Analysis

Global Sensitivity Analysis (GSA) methods are used to quantify the uncertainty in
output of a model w.r.t. the parameters, their individual contributions, or the
contribution of their interactions. The type of GSA method to use depends on
the interest of the user, below we describe the methods available in the suite
at the moment (some more are already in development) and explain what is
the output of each of the methods and what it represents.

## Installation

This functionality does not come standard with DifferentialEquations.jl.
To use this functionality, you must install DiffEqSensitivty.jl:

```julia
]add DiffEqSensitivity
using DiffEqSensitivty
```

## Morris Method

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

## Sobol Method

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

## Regression Method

If a sample of inputs and outputs `` (X^n, Y^n) = 􏰀(X^{i}_1, . . . , X^{i}_d, Y_i)_{i=1..n} ``􏰁
is available, it is possible to fit a linear model explaining the behavior of Y given the
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

## GSA example

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
