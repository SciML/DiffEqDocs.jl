# Random Ordinary Differential Equations

This tutorial will introduce you to the functionality for solving RODEs. Other
introductions can be found by [checking out DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl).

## Example 1: Scalar RODEs

In this example we will solve the equation

```math
du = f(t,u,y)dt
```

where ``f(t,u,W)=2u\sin(W)`` and ``W(t)`` is a Wiener process (Gaussian process).

```julia
using DifferentialEquations
function f(t,u,W)
  2u*sin(W)
end
u0 = 1.00
tspan = (0.0,5.0)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)
```

![intro_rode](../assets/intro_rode.png)

The random process defaults to a Gaussian/Wiener process, so there is nothing
else required here! See the documentation on
[`NoiseProcess`es](../../features/noise_process.html) for details on how to define
other noise proceses.

## Example 2: Systems of RODEs

As with the other problem types, there is an in-place version which is more
efficient for systems. The signature is `f(t,u,W,du)`. For example,

```julia
using DifferentialEquations
function f(t,u,W,du)
  du[1] = 2u[1]*sin(W[1] - W[2])
  du[2] = -2u[2]*cos(W[1] + W[2])
end
u0 = [1.00;1.00]
tspan = (0.0,5.0)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)
```

![rode_system](../assets/rode_system.png)

By default, the size of the noise process matches the size of `u0`. However,
you can use the `rand_prototype` keyword to explicitly set the size of the
random process:

```julia
f = function (t,u,W,du)
  du[1] = -2W[3]*u[1]*sin(W[1] - W[2])
  du[2] = -2u[2]*cos(W[1] + W[2])
end
u0 = [1.00;1.00]
tspan = (0.0,5.0)
prob = RODEProblem(f,u0,tspan,rand_prototype=zeros(3))
sol = solve(prob,RandomEM(),dt=1/100)
```

![noise_choice](../assets/noise_choice.png)
