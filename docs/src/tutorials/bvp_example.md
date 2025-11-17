# Boundary Value Problems

!!! note
    
    This tutorial assumes you have read the [Ordinary Differential Equations tutorial](@ref ode_example).

In this example, we will solve the ODE that satisfies the boundary condition in the form of

```math
\begin{aligned}
\frac{du}{dt} &= f(t, u) \\
g(u) &= \vec{0}
\end{aligned}
```

## Example 1: Simple Pendulum

The concrete example that we are solving is the simple pendulum ``\ddot{u}+\frac{g}{L}sin(u)=0`` on the time interval ``t\in[0,\frac{\pi}{2}]``. First, we need to define the ODE

```@example bvp
import BoundaryValueDiffEq as BVP
import Plots
const g = 9.81
L = 1.0
tspan = (0.0, pi / 2)
function simplependulum!(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
end
```

There are two problem types available:

  - A problem type for general boundary conditions `BVProblem` (including conditions that may be anywhere/ everywhere on the integration interval, aka multi-points BVP).
  - A problem type for boundaries that are specified at the beginning and the end of the integration interval `TwoPointBVProblem`(aka two-points BVP)

The boundary conditions are specified by a function that calculates the residual in-place from the problem solution, such that the residual is $\vec{0}$ when the boundary condition is satisfied.

There are collocation and shooting methods for addressing boundary value problems in DifferentialEquations.jl. We need to use appropriate [available BVP solvers](@ref bvp_solve) to solve `BVProblem`. In this example, we use `MIRK4` to solve the simple pendulum example.

```@example bvp
function bc1!(residual, u, p, t)
    residual[1] = u(pi / 4)[1] + pi / 2 # the solution at the middle of the time span should be -pi/2
    residual[2] = u(pi / 2)[1] - pi / 2 # the solution at the end of the time span should be pi/2
end
bvp1 = BVP.BVProblem(simplependulum!, bc1!, [pi / 2, pi / 2], tspan)
sol1 = BVP.solve(bvp1, BVP.MIRK4(); dt = 0.05)
Plots.plot(sol1)
```

The third argument of `BVProblem` or `TwoPointBVProblem` is the initial guess of the solution, which can be specified as a `Vector`, a `Function` of `t` or solution object from previous solving, in this example the initial guess is set as a `Vector`.

```@example bvp
import OrdinaryDiffEq as ODE
u₀_2 = [-1.6, -1.7] # the initial guess
function bc3!(residual, sol, p, t)
    residual[1] = sol(pi / 4)[1] + pi / 2 # use the interpolation here, since indexing will be wrong for adaptive methods
    residual[2] = sol(pi / 2)[1] - pi / 2
end
bvp3 = BVP.BVProblem(simplependulum!, bc3!, u₀_2, tspan)
sol3 = BVP.solve(bvp3, BVP.Shooting(ODE.Vern7()))
```

The initial guess can also be supplied via a function of `t` or a previous solution type, this is especially handy for parameter analysis.
We changed `u` to `sol` to emphasize the fact that in this case, the boundary condition can be written on the solution object. Thus, all the features on the solution type such as interpolations are available when using both collocation and shooting method. (i.e. you can have a boundary condition saying that the maximum over the interval is `1` using an optimization function on the continuous output).

```@example bvp
Plots.plot(sol3)
```

`TwoPointBVProblem` is operationally the same as `BVProblem` but allows for the solver
to specialize on the common form of being a two-point BVP, i.e. a BVP which only has
boundary conditions at the start and the finish of the time interval.
Defining a similar problem as `TwoPointBVProblem` is shown in the following example:

```@example bvp
function bc2a!(resid_a, u_a, p) # u_a is at the beginning of the time span
    resid_a[1] = u_a[1] + pi / 2 # the solution at the beginning of the time span should be -pi/2
end
function bc2b!(resid_b, u_b, p) # u_b is at the ending of the time span
    resid_b[1] = u_b[1] - pi / 2 # the solution at the end of the time span should be pi/2
end
bvp2 = BVP.TwoPointBVProblem(simplependulum!, (bc2a!, bc2b!), [pi / 2, pi / 2], tspan;
    bcresid_prototype = (zeros(1), zeros(1)))
sol2 = BVP.solve(bvp2, BVP.MIRK4(); dt = 0.05)
Plots.plot(sol2)
```

Note here that `bc2a!` is a boundary condition for the first time point, and `bc2b!` is a boundary condition
for the final time point. `bcresid_prototype` is a prototype array which is passed in order to know the size of
`resid_a` and `resid_b`. In this case, we have one residual term for the start and one for the final time point,
and thus we have `bcresid_prototype = (zeros(1), zeros(1))`.

## Example 2: Directly Solving with Second Order BVP

Suppose we want to solve the second order BVP system which can be formulated as

```math
\begin{cases}
u_1''(x)= u_2(x),\\
\epsilon u_2''(x)=-u_1(x)u_2'(x)- u_3(x)u_3'(x),\\
\epsilon u_3''(x)=u_1'(x)u_3(x)- u_1(x) u_3 '(x)
\end{cases}
```

with initial conditions:

```math
u_1(0) = u_1'(0)= u_1(1)=u_1'(1)=0,u_3(0)=
-1, u_3(1)=1
```

The common way of solving the second order BVP is to define intermediate variables and transform the second order system into first order one, however, DifferentialEquations.jl allows the direct solving of second order BVP system to achieve more efficiency and higher continuity of the numerical solution.

```@example bvp
function f!(ddu, du, u, p, t)
    ϵ = 0.1
    ddu[1] = u[2]
    ddu[2] = (-u[1] * du[2] - u[3] * du[3]) / ϵ
    ddu[3] = (du[1] * u[3] - u[1] * du[3]) / ϵ
end
function bc!(res, du, u, p, t)
    res[1] = u(0.0)[1]
    res[2] = u(1.0)[1]
    res[3] = u(0.0)[3] + 1
    res[4] = u(1.0)[3] - 1
    res[5] = du(0.0)[1]
    res[6] = du(1.0)[1]
end
u0 = [1.0, 1.0, 1.0]
tspan = (0.0, 1.0)
prob = BVP.SecondOrderBVProblem(f!, bc!, u0, tspan)
sol = BVP.solve(
    prob, BVP.MIRKN4(;
        jac_alg = BVP.BVPJacobianAlgorithm(BVP.AutoForwardDiff())); dt = 0.01)
```

## Example 3: Semi-Explicit Boundary Value Differential-Algebraic Equations

Consider a semi-explicit boundary value differential-algebraic equation formulated as

```math
\begin{cases}
x_1'=(\epsilon+x_2-p_2(t))y+p_1'(t) \\
x_2'=p_2'(t) \\
x_3'=y \\
0=(x_1-p_1(t))(y-e^t)
\end{cases}
```

with boundary conditions

```math
x_1(0)=0,x_3(0)=1,x_2(1)=\sin(1)
```

We need to choose the Ascher methods for solving BVDAEs.

```@example bvp
function f!(du, u, p, t)
    e = 2.7
    du[1] = (1 + u[2] - sin(t)) * u[4] + cos(t)
    du[2] = cos(t)
    du[3] = u[4]
    du[4] = (u[1] - sin(t)) * (u[4] - e^t)
end
function bc!(res, u, p, t)
    res[1] = u[1]
    res[2] = u[3] - 1
    res[3] = u[2] - sin(1.0)
end
u0 = [0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 1.0)
fun = BVP.BVPFunction(f!, bc!, mass_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0])
prob = BVP.BVProblem(fun, u0, tspan)
sol = BVP.solve(prob,
    BVP.Ascher4(; zeta = [0.0, 0.0, 1.0], jac_alg = BVP.BVPJacobianAlgorithm(BVP.AutoForwardDiff()));
    dt = 0.01)
```
