# Differential Algebraic Equations

This tutorial will introduce you to the functionality for solving DAEs. Other
introductions can be found by [checking out DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl). This tutorial assumes you have read the [Ordinary Differential Equations tutorial](ode_example.html).

In this example we will solve the implicit ODE equation

```math
f(du,u,p,t) = 0
```

where `f` is the a variant of the Roberts equation. This equation is of
the form

```math
\begin{align}
du &= f(u,p,t) \\
 0 &= g(u,p,t) \\
 \end{align}
```

or is also known as a constrained differential equation where `g` is the constraint
equation. The Robertson model can be written in the form:

```math
\begin{align}
dy_1 &= -0.04y₁ + 10^4 y_2 y_3 \\
dy_2 &= 0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
1 &=  y_{1} + y_{2} + y_{3} \\
\end{align}
```

with initial conditions ``y_1(0) = 1``, ``y_2(0) = 0``, ``y_3(0) = 0``,
``dy_1 = - 0.04``, ``dy_2 = 0.04``, and ``dy_3 = 0.0``.

The workflow for DAEs is the same as for the other types of equations, where all
you need to know is how to define the problem. A DAEProblem is specified by defining
an in-place update `f(out,du,u,p,t)` which uses the values to mutate `out` as the
output. To makes this into a DAE, we move all of the variables to one side.
Thus we can define the function:

```julia
function f(out,du,u,p,t)
  out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
  out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
  out[3] = u[1] + u[2] + u[3] - 1.0
end
```

with initial conditions

```julia
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0,100000.0)
```

and make the DAEProblem:

```julia
using DifferentialEquations
differential_vars = [true,true,false]
prob = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)
```

`differential_vars` is an option which states which of the variables are differential,
i.e. not purely algebraic (which means that their derivative shows up in the residual
equations). This is required for the algorithm to be able to find consistant initial
conditions. Notice that the first two variables are determined by their changes, but
the last is simply determined by the conservation equation. Thus we use
`differential_vars = [true,true,false]`.

As with the other DifferentialEquations problems, the commands are then to solve
and plot. Here we will use the IDA solver from Sundials:

```julia
using Sundials
sol = solve(prob,IDA())
```

In order to clearly see all the features of this solution, it should be plotted
on a logarithmic scale. We'll also plot each on a different subplot to allow
scaling the y-axis appropriately.

```julia
using Plots; plotly() # Using the Plotly backend
plot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))
```

This gives the following plot:

![IntroDAEPlot](../assets/intro_dae_plot.png)
