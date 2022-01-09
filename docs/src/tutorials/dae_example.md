# Differential Algebraic Equations

Differential Algebraic Equations (DAEs) are differential equations which have
constraint equations on their evolution. This tutorial will introduce you to
the functionality for solving differential algebraic equations (DAEs). Other
introductions can be found by
[checking out SciMLTutorials.jl](https://github.com/JuliaDiffEq/SciMLTutorials.jl).

!!! note

    This tutorial assumes you have read the [Ordinary Differential Equations tutorial](@ref ode_example).

## Mass-Matrix Differential-Algebraic Equations (DAEs)

Instead of just defining an ODE as ``u' = f(u,p,t)``, it can be common to express
the differential equation in the form with a mass matrix:

```math
Mu' = f(u,p,t)
```

where ``M`` is known as the mass matrix. Let's solve the Robertson equation.
In previous tutorials we wrote this equation as:

```math
\begin{aligned}
dy_1 &= -0.04 y_1 + 10^4 y_2 y_3 \\
dy_2 &=  0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
dy_3 &= 3*10^7 y_{2}^2 \\
\end{aligned}
```

But we can instead write this with a conservation relation:

```math
\begin{aligned}
\frac{dy_1}{dt} &= -0.04 y_1 + 10^4 y_2 y_3 \\
\frac{dy_2}{dt} &=  0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
1 &=  y_{1} + y_{2} + y_{3} \\
\end{aligned}
```

In this form, we can write this as a mass matrix ODE where ``M`` is singular
(this is another form of a differential-algebraic equation (DAE)). Here, the
last row of `M` is just zero. We can implement this form as:

```julia
using DifferentialEquations
function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁ + k₃*y₂*y₃
  du[2] =  k₁*y₁ - k₃*y₂*y₃ - k₂*y₂^2
  du[3] =  y₁ + y₂ + y₃ - 1
  nothing
end
M = [1. 0  0
     0  1. 0
     0  0  0]
f = ODEFunction(rober,mass_matrix=M)
prob_mm = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol = solve(prob_mm,Rodas5(),reltol=1e-8,abstol=1e-8)

plot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))
```

![IntroDAEPlot](../assets/intro_dae_plot.png)

!!! note
    If your mass matrix is singular, i.e. your system is a DAE, then you
    need to make sure you choose
    [a solver that is compatible with DAEs](@ref dae_solve_full)

## Implicitly-Defined Differential-Algebraic Equations (DAEs)

In this example we will solve the Robertson equation in its implicit form:

```math
f(du,u,p,t) = 0
```

This equation is a DAE of the form:

```math
\begin{aligned}
\frac{du}{dt} &= f(u,p,t) \\
 0 &= g(u,p,t) \\
 \end{aligned}
```

which is also known as a constrained differential equation, where `g` is the constraint
equation. The Robertson model can be written in the form:

```math
\begin{aligned}
\frac{dy_1}{dt} &= -0.04y₁ + 10^4 y_2 y_3 \\
\frac{dy_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3*10^7 y_{2}^2 \\
1 &=  y_{1} + y_{2} + y_{3} \\
\end{aligned}
```

with initial conditions ``y_1(0) = 1``, ``y_2(0) = 0``, ``y_3(0) = 0``,
``dy_1 = - 0.04``, ``dy_2 = 0.04``, and ``dy_3 = 0.0``.

The workflow for DAEs is the same as for the other types of equations, where all
you need to know is how to define the problem. A `DAEProblem` is specified by defining
an in-place update `f(out,du,u,p,t)` which uses the values to mutate `out` as the
output. To makes this into a DAE, we move all of the variables to one side.
Thus, we can define the function:

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

and make the `DAEProblem`:

```julia
using DifferentialEquations
differential_vars = [true,true,false]
prob = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)
```

`differential_vars` is an option which states which of the variables are differential,
i.e. not purely algebraic (which means that their derivative shows up in the residual
equations). This is required for the algorithm to be able to find consistent initial
conditions. Notice that the first two variables are determined by their changes, but
the last is simply determined by the conservation equation. Thus, we use
`differential_vars = [true,true,false]`.

As with the other DifferentialEquations problems, the commands are then to solve
and plot. Here we will use the IDA solver from Sundials:

```julia
using Sundials
sol = solve(prob,IDA())
```

In order to clearly see all the features of this solution, it should be plotted
on a logarithmic scale. We'll also plot each on a different subplot, to allow
scaling the y-axis appropriately.

```julia
using Plots
plot(sol, xscale=:log10, tspan=(1e-6, 1e5), layout=(3,1))
```

This gives the following plot:

![IntroDAEPlot](../assets/intro_dae_plot.png)
