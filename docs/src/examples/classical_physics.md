# Classical Physics Models

If you're getting some cold feet to jump in to DiffEq land, here are some handcrafted differential equations mini problems to hold your hand along the beginning of your journey.

## First order linear ODE

#### Radioactive Decay of Carbon-14

$$f(t,u) = \frac{du}{dt}$$

The Radioactive decay problem is the first order linear ODE problem of an exponential with a negative coefficient, which represents the half-life of the process in question. Should the coefficient be positive, this would represent a population growth equation.

```@example physics
import OrdinaryDiffEq as ODE, Plots
Plots.gr()

#Half-life of Carbon-14 is 5,730 years.
t½ = 5.730

#Setup
u₀ = 1.0
tspan = (0.0, 30.0)

#Define the problem
radioactivedecay(u, p, t) = -log(2) / t½ * u

#Pass to solver
prob = ODE.ODEProblem(radioactivedecay, u₀, tspan)
sol = ODE.solve(prob, ODE.Tsit5())

#Plot
Plots.plot(sol, linewidth = 2, title = "Carbon-14 half-life",
    xaxis = "Time in thousands of years", yaxis = "Ratio left",
    label = "Numerical Solution")
Plots.plot!(sol.t, t -> 2^(-t / t½), lw = 3, ls = :dash, label = "Analytical Solution")
```

## Second Order Linear ODE

#### Simple Harmonic Oscillator

Another classical example is the harmonic oscillator, given by:

$$\ddot{x} + \omega^2 x = 0$$

with the known analytical solution

$$\begin{align*}
x(t) &= A\cos(\omega t - \phi) \\
v(t) &= -A\omega\sin(\omega t - \phi),
\end{align*}$$

where

$$A = \sqrt{c_1 + c_2} \qquad\text{and}\qquad \tan \phi = \frac{c_2}{c_1}$$

with $c_1, c_2$ constants determined by the initial conditions such that
$c_1$ is the initial position and $\omega c_2$ is the initial velocity.

Instead of transforming this to a system of ODEs to solve with `ODEProblem`,
we can use `SecondOrderODEProblem` as follows.

```@example physics
# Simple Harmonic Oscillator Problem
import OrdinaryDiffEq as ODE, Plots

#Parameters
ω = 1

#Initial Conditions
x₀ = [0.0]
dx₀ = [π / 2]
tspan = (0.0, 2π)

ϕ = atan((dx₀[1] / ω) / x₀[1])
A = √(x₀[1]^2 + dx₀[1]^2)

#Define the problem
function harmonicoscillator(ddu, du, u, ω, t)
    ddu .= -ω^2 * u
end

#Pass to solvers
prob = ODE.SecondOrderODEProblem(harmonicoscillator, dx₀, x₀, tspan, ω)
sol = ODE.solve(prob, ODE.DPRKN6())

#Plot
Plots.plot(sol, idxs = [2, 1], linewidth = 2, title = "Simple Harmonic Oscillator",
    xaxis = "Time", yaxis = "Elongation", label = ["x" "dx"])
Plots.plot!(t -> A * cos(ω * t - ϕ), lw = 3, ls = :dash, label = "Analytical Solution x")
Plots.plot!(t -> -A * ω * sin(ω * t - ϕ), lw = 3, ls = :dash, label = "Analytical Solution dx")
```

Note that the order of the variables (and initial conditions) is `dx`, `x`.
Thus, if we want the first series to be `x`, we have to flip the order with `vars=[2,1]`.

## Second Order Non-linear ODE

#### Simple Pendulum

We will start by solving the pendulum problem. In the physics class, we often solve this problem by small angle approximation, i.e. $ sin(\theta) \approx \theta$, because otherwise, we get an elliptic integral which doesn't have an analytic solution. The linearized form is

$$\ddot{\theta} + \frac{g}{L}{\theta} = 0$$

But we have numerical ODE solvers! Why not solve the *real* pendulum?

$$\ddot{\theta} + \frac{g}{L}{\sin(\theta)} = 0$$

Notice that now we have a second order ODE.
In order to use the same method as above, we need to transform it into a system
of first order ODEs by employing the notation $d\theta = \dot{\theta}$.

$$\begin{align*}
&\dot{\theta} = d{\theta} \\
&\dot{d\theta} = - \frac{g}{L}{\sin(\theta)}
\end{align*}$$

```@example physics
# Simple Pendulum Problem
import OrdinaryDiffEq as ODE, Plots

#Constants
const g = 9.81
L = 1.0

#Initial Conditions
u₀ = [0, π / 2]
tspan = (0.0, 6.3)

#Define the problem
function simplependulum(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
end

#Pass to solvers
prob = ODE.ODEProblem(simplependulum, u₀, tspan)
sol = ODE.solve(prob, ODE.Tsit5())

#Plot
Plots.plot(sol, linewidth = 2, title = "Simple Pendulum Problem", xaxis = "Time",
    yaxis = "Height", label = ["\\theta" "d\\theta"])
```

So now we know that behaviour of the position versus time. However, it will be useful to us to look at the phase space of the pendulum, i.e., and representation of all possible states of the system in question (the pendulum) by looking at its velocity and position. Phase space analysis is ubiquitous in the analysis of dynamical systems, and thus we will provide a few facilities for it.

```@example physics
p = Plots.plot(sol, vars = (1, 2), xlims = (-9, 9), title = "Phase Space Plot",
    xaxis = "Angular position", yaxis = "Angular velocity", leg = false)
function phase_plot(prob, u0, p, tspan = 2pi)
    _prob = ODE.ODEProblem(prob.f, u0, (0.0, tspan))
    sol = ODE.solve(_prob, ODE.Vern9()) # Use Vern9 solver for higher accuracy
    Plots.plot!(p, sol, idxs = (1, 2), xlims = nothing, ylims = nothing)
end
for i in (-4pi):(pi / 2):(4π)
    for j in (-4pi):(pi / 2):(4π)
        phase_plot(prob, [j, i], p)
    end
end
Plots.plot(p, xlims = (-9, 9))
```

#### Double Pendulum

A more complicated example is given by the double pendulum. The equations governing
its motion are given by the following (taken from this [Stack Overflow question](https://mathematica.stackexchange.com/questions/40122/help-to-plot-poincar%C3%A9-section-for-double-pendulum))

$$\frac{d}{dt}
\begin{pmatrix}
\alpha \\ l_\alpha \\ \beta \\ l_\beta
\end{pmatrix}=
\begin{pmatrix}
2\frac{l_\alpha - (1+\cos\beta)l_\beta}{3-\cos 2\beta} \\
-2\sin\alpha - \sin(\alpha + \beta) \\
2\frac{-(1+\cos\beta)l_\alpha + (3+2\cos\beta)l_\beta}{3-\cos2\beta}\\
-\sin(\alpha+\beta) - 2\sin(\beta)\frac{(l_\alpha-l_\beta)l_\beta}{3-\cos2\beta} + 2\sin(2\beta)\frac{l_\alpha^2-2(1+\cos\beta)l_\alpha l_\beta + (3+2\cos\beta)l_\beta^2}{(3-\cos2\beta)^2}
\end{pmatrix}$$

```@example physics
#Double Pendulum Problem
import OrdinaryDiffEq as ODE, Plots

#Constants and setup
const m₁, m₂, L₁, L₂ = 1, 2, 1, 2
initial = [0, π / 3, 0, 3pi / 5]
tspan = (0.0, 50.0)

#Convenience function for transforming from polar to Cartesian coordinates
function polar2cart(sol; dt = 0.02, l1 = L₁, l2 = L₂, vars = (2, 4))
    u = sol.t[1]:dt:sol.t[end]

    p1 = l1 * map(x -> x[vars[1]], sol.(u))
    p2 = l2 * map(y -> y[vars[2]], sol.(u))

    x1 = l1 * sin.(p1)
    y1 = l1 * -cos.(p1)
    (u, (x1 + l2 * sin.(p2),
        y1 - l2 * cos.(p2)))
end

#Define the Problem
function double_pendulum(xdot, x, p, t)
    xdot[1] = x[2]
    xdot[2] = -((g * (2 * m₁ + m₂) * sin(x[1]) +
                 m₂ * (g * sin(x[1] - 2 * x[3]) +
                  2 * (L₂ * x[4]^2 + L₁ * x[2]^2 * cos(x[1] - x[3])) * sin(x[1] - x[3]))) /
                (2 * L₁ * (m₁ + m₂ - m₂ * cos(x[1] - x[3])^2)))
    xdot[3] = x[4]
    xdot[4] = (((m₁ + m₂) * (L₁ * x[2]^2 + g * cos(x[1])) +
                L₂ * m₂ * x[4]^2 * cos(x[1] - x[3])) * sin(x[1] - x[3])) /
              (L₂ * (m₁ + m₂ - m₂ * cos(x[1] - x[3])^2))
end

#Pass to Solvers
double_pendulum_problem = ODE.ODEProblem(double_pendulum, initial, tspan)
sol = ODE.solve(double_pendulum_problem, ODE.Vern7(), abstol = 1e-10, dt = 0.05);
```

```@example physics
#Obtain coordinates in Cartesian Geometry
ts, ps = polar2cart(sol, l1 = L₁, l2 = L₂, dt = 0.01)
Plots.plot(ps...)
```

##### Poincaré section

In this case, the phase space is 4 dimensional, and it cannot be easily visualized.
Instead of looking at the full phase space, we can look at Poincaré sections,
which are sections through a higher-dimensional phase space diagram.
This helps to understand the dynamics of interactions and is wonderfully pretty.

The Poincaré section in this is given by the collection of $(β,l_β)$ when $α=0$ and $\frac{dα}{dt}>0$.

```@example physics
#Constants and setup
import OrdinaryDiffEq as ODE
initial2 = [0.01, 0.005, 0.01, 0.01]
tspan2 = (0.0, 500.0)

#Define the problem
function double_pendulum_hamiltonian(udot, u, p, t)
    α = u[1]
    lα = u[2]
    β = u[3]
    lβ = u[4]
    udot .= [2(lα - (1 + cos(β))lβ) / (3 - cos(2β)),
        -2sin(α) - sin(α + β),
        2(-(1 + cos(β))lα + (3 + 2cos(β))lβ) / (3 - cos(2β)),
        -sin(α + β) - 2sin(β) * (((lα - lβ)lβ) / (3 - cos(2β))) +
        2sin(2β) * ((lα^2 - 2(1 + cos(β))lα * lβ + (3 + 2cos(β))lβ^2) / (3 - cos(2β))^2)]
end

# Construct a ContiunousCallback
condition(u, t, integrator) = u[1]
affect!(integrator) = nothing
cb = ODE.ContinuousCallback(condition, affect!, nothing,
    save_positions = (true, false))

# Construct Problem
poincare = ODE.ODEProblem(double_pendulum_hamiltonian, initial2, tspan2)
sol2 = ODE.solve(poincare, ODE.Vern9(), save_everystep = false, save_start = false,
    save_end = false, callback = cb, abstol = 1e-16, reltol = 1e-16)

function poincare_map(prob, u₀, p; callback = cb)
    _prob = ODE.ODEProblem(prob.f, u₀, prob.tspan)
    sol = ODE.solve(_prob, ODE.Vern9(), save_everystep = false, save_start = false,
        save_end = false, callback = cb, abstol = 1e-16, reltol = 1e-16)
    Plots.scatter!(p, sol, idxs = (3, 4), markersize = 3, msw = 0)
end
```

```@example physics
lβrange = -0.02:0.0025:0.02
p = Plots.scatter(sol2, idxs = (3, 4), leg = false, markersize = 3, msw = 0)
for lβ in lβrange
    poincare_map(poincare, [0.01, 0.01, 0.01, lβ], p)
end
Plots.plot(p, xlabel = "\\beta", ylabel = "l_\\beta", ylims = (0, 0.03))
```

#### Hénon-Heiles System

The Hénon-Heiles potential occurs when non-linear motion of a star around a galactic center, with the motion restricted to a plane.

$$\begin{align}
\frac{d^2x}{dt^2}&=-\frac{\partial V}{\partial x}\\
\frac{d^2y}{dt^2}&=-\frac{\partial V}{\partial y}
\end{align}$$

where

$$V(x,y)={\frac {1}{2}}(x^{2}+y^{2})+\lambda \left(x^{2}y-{\frac {y^{3}}{3}}\right).$$

We pick $\lambda=1$ in this case, so

$$V(x,y) = \frac{1}{2}(x^2+y^2+2x^2y-\frac{2}{3}y^3).$$

Then the total energy of the system can be expressed by

$$E = T+V = V(x,y)+\frac{1}{2}(\dot{x}^2+\dot{y}^2).$$

The total energy should conserve as this system evolves.

```@example physics
import OrdinaryDiffEq as ODE, Plots

#Setup
initial = [0.0, 0.1, 0.5, 0]
tspan = (0, 100.0)

#Remember, V is the potential of the system and T is the Total Kinetic Energy, thus E will
#the total energy of the system.
V(x, y) = 1 // 2 * (x^2 + y^2 + 2x^2 * y - 2 // 3 * y^3)
E(x, y, dx, dy) = V(x, y) + 1 // 2 * (dx^2 + dy^2);

#Define the function
function Hénon_Heiles(du, u, p, t)
    x = u[1]
    y = u[2]
    dx = u[3]
    dy = u[4]
    du[1] = dx
    du[2] = dy
    du[3] = -x - 2x * y
    du[4] = y^2 - y - x^2
end

#Pass to solvers
prob = ODE.ODEProblem(Hénon_Heiles, initial, tspan)
sol = ODE.solve(prob, ODE.Vern9(), abstol = 1e-16, reltol = 1e-16);
```

```@example physics
# Plot the orbit
Plots.plot(sol, idxs = (1, 2), title = "The orbit of the Hénon-Heiles system", xaxis = "x",
    yaxis = "y", leg = false)
```

```@example physics
#Optional Sanity check - what do you think this returns and why?
@show sol.retcode

#Plot -
Plots.plot(sol, idxs = (1, 3), title = "Phase space for the Hénon-Heiles system",
    xaxis = "Position", yaxis = "Velocity")
Plots.plot!(sol, idxs = (2, 4), leg = false)
```

```@example physics
#We map the Total energies during the time intervals of the solution (sol.u here) to a new vector
#pass it to the plotter a bit more conveniently
energy = map(x -> E(x...), sol.u)

#We use @show here to easily spot erratic behavior in our system by seeing if the loss in energy was too great.
@show ΔE = energy[1] - energy[end]

#Plot
Plots.plot(sol.t, energy .- energy[1], title = "Change in Energy over Time",
    xaxis = "Time in iterations", yaxis = "Change in Energy")
```

##### Symplectic Integration

To prevent energy drift, we can instead use a symplectic integrator. We can directly define and solve the `SecondOrderODEProblem`:

```@example physics
function HH_acceleration!(dv, v, u, p, t)
    x, y = u
    dx, dy = dv
    dv[1] = -x - 2x * y
    dv[2] = y^2 - y - x^2
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
prob = ODE.SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol2 = ODE.solve(prob, ODE.KahanLi8(), dt = 1 / 10);
```

Notice that we get the same results:

```@example physics
# Plot the orbit
Plots.plot(sol2, idxs = (3, 4), title = "The orbit of the Hénon-Heiles system", xaxis = "x",
    yaxis = "y", leg = false)
```

```@example physics
Plots.plot(sol2, idxs = (3, 1), title = "Phase space for the Hénon-Heiles system",
    xaxis = "Position", yaxis = "Velocity")
Plots.plot!(sol2, idxs = (4, 2), leg = false)
```

but now the energy change is essentially zero:

```@example physics
energy = map(x -> E(x[3], x[4], x[1], x[2]), sol2.u)
#We use @show here to easily spot erratic behaviour in our system by seeing if the loss in energy was too great.
@show ΔE = energy[1] - energy[end]

#Plot
Plots.plot(sol2.t, energy .- energy[1], title = "Change in Energy over Time",
    xaxis = "Time in iterations", yaxis = "Change in Energy")
```

And let's try to use a Runge-Kutta-Nyström solver to solve this. Note that Runge-Kutta-Nyström isn't symplectic.

```@example physics
sol3 = ODE.solve(prob, ODE.DPRKN6());
energy = map(x -> E(x[3], x[4], x[1], x[2]), sol3.u)
@show ΔE = energy[1] - energy[end]
Plots.gr()
Plots.plot(sol3.t, energy .- energy[1], title = "Change in Energy over Time",
    xaxis = "Time in iterations", yaxis = "Change in Energy")
```

Note that we are using the `DPRKN6` solver at `reltol=1e-3` (the default), yet it has a smaller energy variation than `Vern9` at `abstol=1e-16, reltol=1e-16`. Therefore, using specialized solvers to solve its particular problem is very efficient.
