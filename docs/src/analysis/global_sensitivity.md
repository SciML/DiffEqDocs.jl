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
using DiffEqSensitivity
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
effects for each parameter which takes into consideration all the interactions and its
individual contribution.

`morris_effects = gsa(f,Morris(),param_range)`


Here, `f` is just the model (as a julia function 
`f(input_vector) -> output_vector` or a `DEProblem`) you want to
run the analysis on. The `Morris` object signifies the method to be used and contains
the following fields that can be passed by user to adjust the parameter sampling:

  1. `p_steps` - Decides the value of ``\Delta`` in the elementary effects calculation.
  2. `relative_scale` - The elementary effects are calculated with the assumption that
the parameters lie in the range `[0,1]` but as this is not always the case
scaling is used to get more informative, scaled effects.
  3. `total_num_trajectory`, `num_trajectory` - The total number of design matrices that are 
  generated out of which `num_trajectory` matrices with the highest spread are used in calculation.
  4. `len_design_mat` - The size of a design matrix.

  `param_range` requires an array of 2-tuples with the lower bound
and the upper bound.

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

`sobol_indices = gsa(f,Sobol(),A,B;batch=false,Ei_estimator = :Jansen1999)`

The `Sobol` object has as its fields the `order` of the indices to be estimated. 
The `Ei_estimator` kwarg can take `:Homma1996`, `:Sobol2007` and `:Jansen1999` 
which signify the different ways of estimating the indices.   

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
function f(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = -3*u[2] + u[1]*u[2]
end
u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5,1.0]
prob = ODEProblem(f,u0,tspan,p)
t = collect(range(0, stop=10, length=200))
```
For Morris Method

```julia
f1 = function (p)
        prob1 = remake(prob;p=p)
        Array(solve(prob1,Tsit5();saveat=t))
    end
m = gsa(f1,Morris(total_num_trajectory=1000,num_trajectory=150),[[1,5],[0.5,5]])
```
Let's get the means and variances from the `MorrisResult` struct.

```julia
m.means
2-element Array{Array{Float64,2},1}:
 [0.0 0.051096 … 12.5698 12.401; 0.0 0.00115346 … 5.34451 5.63429]
 [0.0 0.048616 … 2.31299 2.39601; 0.0 0.00111577 … 2.5087 2.4553] 

m.variances
2-element Array{Array{Float64,2},1}:
 [0.0 1.92803e-5 … 160.412 157.503; 0.0 4.77857e-9 … 84.1636 88.2436]
 [0.0 1.70976e-5 … 21.9696 27.3176; 0.0 4.3237e-9 … 86.6352 77.7245] 
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
f1 = function (p)
        prob1 = remake(prob;p=p)
        Array(solve(prob1,Tsit5();saveat=t))[end]
    end

N = 10000
lb = [1.0, 0.5]
ub = [5.0, 5.0]
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(N,lb,ub,sampler)
sobol_result = gsa(f1,Sobol(),A,B)
```

We plot the first order and total order Sobol Indices for some timepoints for each of the parameters (`a` and `b`).

```julia

p1 = bar(["a","b"],[sobol_result.ST[1],sobol_result.ST[2]],color=[:red,:blue],title="Total Order Indices",legend=false)
p2 = bar(["a","b"],[sobol_result.S1[1],sobol_result.S1[2]],color=[:red,:blue],title="First Order Indices",legend=false)
plo = plot(p1,p2,layout=(2,1),size=(600,500))

```
![sobolplot](../assets/sobolbars.png)

Here we plot the Sobol indices of first order and the total Sobol indices for the parameters `a` and `b`. The plots are obtained by getting the Sobol Indices at the final time point of the second dependent variable `y(t)` from the 200-length sensitivities over the entire time span. The length of the bar represents the quantification of the sensitivity of the output to that parameter and hence for the final time point you can say that `y(t)` is more sensitive to `b`, also you can observe how the relative difference between `a` and `b` is larger in the first order than the total order indices, this tells us that most of the contribution of `a` to `y(t)` arises from interactions and it's individual non-interaction contribution is significantly lesser than `b` and vice-versa for `b` as it's first order plot indicates quite high value.
