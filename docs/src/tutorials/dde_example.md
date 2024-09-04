# Delay Differential Equations

This tutorial will introduce you to the functionality for solving delay differential
equations.

!!! note
    
    This tutorial assumes you have read the [Ordinary Differential Equations tutorial](@ref ode_example).

Delay differential equations are equations which have a delayed argument. To allow
for specifying the delayed argument, the function definition for a delay differential
equation is expanded to include a history function `h(p, t)` which uses interpolations
throughout the solution's history to form a continuous extension of the solver's
past and depends on parameters `p` and time `t`. The function signature for a delay
differential equation is `f(u, h, p, t)` for not in-place computations, and
`f(du, u, h, p, t)` for in-place computations.

In this example, we will solve [a model of breast cancer growth kinetics](https://idp.nature.com/authorize?response_type=cookie&client_id=grover&redirect_uri=https%3A%2F%2Fwww.nature.com%2Farticles%2Fsrep02473):

```math
\begin{aligned}
dx_{0} &= \frac{v_{0}}{1+\beta_{0}\left(x_{2}(t-\tau)\right)^{2}}\left(p_{0}-q_{0}\right)x_{0}(t)-d_{0}x_{0}(t)\\
dx_{1} &= \frac{v_{0}}{1+\beta_{0}\left(x_{2}(t-\tau)\right)^{2}}\left(1-p_{0}+q_{0}\right)x_{0}(t)\\
       &+ \frac{v_{1}}{1+\beta_{1}\left(x_{2}(t-\tau)\right)^{2}}\left(p_{1}-q_{1}\right)x_{1}(t)-d_{1}x_{1}(t)\\
dx_{2} &= \frac{v_{1}}{1+\beta_{1}\left(x_{2}(t-\tau)\right)^{2}}\left(1-p_{1}+q_{1}\right)x_{1}(t)-d_{2}x_{2}(t)
\end{aligned}
```

For this problem, we note that ``\tau`` is constant, and thus we can use a method
which exploits this behavior. We first write out the equation using the appropriate
function signature. Most of the equation writing is the same, though we use the
history function by first interpolating and then choosing the components. Thus,
the `i`th component at time `t-tau` is given by `h(p, t-tau)[i]`. Components with
no delays are written as in the ODE.

Thus, the function for this model is given by:

```@example dde
using DifferentialEquations
function bc_model(du, u, h, p, t)
    p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau = p
    hist3 = h(p, t - tau)[3]
    du[1] = (v0 / (1 + beta0 * (hist3^2))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (hist3^2))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (hist3^2))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (hist3^2))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end
```

Now we build a `DDEProblem`. The signature

```julia
prob = DDEProblem(f, u0, h, tspan, p = SciMLBase.NullParameters();
    constant_lags = [], dependent_lags = [], kwargs...)
```

is very similar to ODEs, where we now have to give the lags and a function `h`.
`h` is the history function that declares what the values were before the time
the model starts. Here we will assume that for all time before `t0` the values were 1
and define `h` as an out-of-place function:

```@example dde
h(p, t) = ones(3)
```

To use the constant lag model, we have to declare the lags. Here we will use `tau=1`.

```@example dde
tau = 1
lags = [tau]
```

Next, we choose to solve on the timespan `(0.0,10.0)` and create the problem type:

```@example dde
p0 = 0.2;
q0 = 0.3;
v0 = 1;
d0 = 5;
p1 = 0.2;
q1 = 0.3;
v1 = 1;
d1 = 1;
d2 = 1;
beta0 = 1;
beta1 = 1;
p = (p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau)
tspan = (0.0, 10.0)
u0 = [1.0, 1.0, 1.0]

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)
```

An efficient way to solve this problem (given the constant lags) is with the
MethodOfSteps solver. Through the magic that is Julia, it translates an OrdinaryDiffEq.jl
ODE solver method into a method for delay differential equations, which is highly
efficient due to sweet compiler magic. A good choice is the order 5 method `Tsit5()`:

```@example dde
alg = MethodOfSteps(Tsit5())
```

For lower tolerance solving, one can use the `BS3()` algorithm to good
effect (this combination is similar to the MATLAB `dde23`, but more efficient
tableau), and for high tolerances the `Vern6()` algorithm will give a 6th order
solution.

To solve the problem with this algorithm, we do the same thing we'd do with other
methods on the common interface:

```@example dde
sol = solve(prob, alg)
```

Note that everything available to OrdinaryDiffEq.jl can be used here, including
event handling and other callbacks. The solution object has the same interface
as for ODEs. For example, we can use the same plot recipes to view the results:

```@example dde
using Plots
plot(sol)
```

#### Speeding Up Interpolations with Idxs

We can speed up the previous problem in two different ways. First of all, if we
need to interpolate multiple values from a previous time, we can use the in-place
form for the history function `h(out, p, t)` which writes the output to `out`. In this
case, we must supply the history initial conditions as in-place as well. For the
previous example, that's simply

```@example dde
h(out, p, t) = (out .= 1.0)
```

and then our DDE is:

```@example dde
const out = zeros(3) # Define a cache variable
function bc_model(du, u, h, p, t)
    h(out, p, t - tau) # updates out to be the correct history function
    du[1] = (v0 / (1 + beta0 * (out[3]^2))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (out[3]^2))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (out[3]^2))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (out[3]^2))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end
```

However, we can do something even slicker in most cases. We only ever needed to
interpolate past values at index 3. Instead of generating a bunch of arrays,
we can instead ask specifically for that value by passing the keyword `idxs = 3`.
The DDE function `bc_model` is now:

```@example dde
function bc_model(du, u, h, p, t)
    u3_past_sq = h(p, t - tau; idxs = 3)^2
    du[1] = (v0 / (1 + beta0 * (u3_past_sq))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (u3_past_sq))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (u3_past_sq))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (u3_past_sq))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end
```

Note that this requires that we define the historical values

```@example dde
h(p, t; idxs = nothing) = typeof(idxs) <: Number ? 1.0 : ones(3)
```

where `idxs` can be an integer for which variable in the history to compute,
and here for any number `idxs` we give back `1.0`. Note that if we wanted to use
past values of the `i`th derivative, then we would call the history function
`h(p, t, Val{i})` in our DDE function and would have to define a dispatch like

```@example dde
h(p, t, ::Type{Val{1}}) = zeros(3)
```

to say that derivatives before `t0` are zero for any index. Again, we could
use an in-place function instead or only compute specific indices by passing
an `idxs` keyword.

The functional forms for the history function are also discussed
[on the DDEProblem page](@ref dde_prob).

### Undeclared Delays and State-Dependent Delays via Residual Control

You might have noticed DifferentialEquations.jl allows you to solve problems
with undeclared delays, since you can interpolate `h` at any value. This is
a feature, but use it with caution. Undeclared delays can increase the error
in the solution. It's recommended that you use a method with a residual control,
such as `MethodOfSteps(RK4())` whenever there are undeclared delays. With this,
you can use interpolated derivatives, solve functional differential equations
by using quadrature on the interpolant, etc. However, note that residual control
solves with a low level of accuracy, so the tolerances should be made very small,
and the solution should not be trusted for more than 2-3 decimal places.

Note: `MethodOfSteps(RK4())` with undeclared delays is similar to MATLAB's
`ddesd`. Thus, for example, the following is similar to solving the example
from above with residual control:

```@example dde
prob = DDEProblem(bc_model, u0, h, tspan)
alg = MethodOfSteps(RK4())
sol = solve(prob, alg)
```

Note that this method can solve problems with state-dependent delays.

### State-Dependent Delay Discontinuity Tracking

State-dependent delays are problems where the delay is allowed to be a function
of the current state. They can be more efficiently solved with discontinuity
tracking. To do this, in DifferentialEquations.jl, requires passing lag functions
`g(u,p,t)` as keyword `dependent_lags` to the `DDEProblem` definition. Other than
that, everything else is the same, and one solves that problem using the common
interface.

We can solve the above problem with dependent delay tracking by declaring the
dependent lags and solving with a `MethodOfSteps` algorithm:

```@example dde
prob = DDEProblem(bc_model, u0, h, tspan; dependent_lags = ((u, p, t) -> tau,))
alg = MethodOfSteps(Tsit5())
sol = solve(prob, alg)
```

Here, we treated the single lag `t-tau` as a state-dependent delay. Of course, you
can then replace that tuple of functions with whatever functions match your lags.
