# Finding Maxima and Minima of ODEs Solutions

### Setup

In this tutorial, we will show how to use
[Optimization.jl](https://docs.sciml.ai/Optimization/stable/) to find the maxima and minima
of solutions. Let's take a look at the double pendulum:

```@example minmax
#Constants and setup
import OrdinaryDiffEq as ODE
initial = [0.01, 0.01, 0.01, 0.01]
tspan = (0.0, 100.0)

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

#Pass to solvers
poincare = ODE.ODEProblem(double_pendulum_hamiltonian, initial, tspan)
```

```@example minmax
sol = ODE.solve(poincare, ODE.Tsit5())
```

In time, the solution looks like:

```@example minmax
import Plots;
Plots.gr();
Plots.plot(sol, vars = [(0, 3), (0, 4)], leg = false, plotdensity = 10000)
```

while it has the well-known phase-space plot:

```@example minmax
Plots.plot(sol, vars = (3, 4), leg = false)
```

### Local Optimization

Let's find out what some of the local maxima and minima are. Optim.jl can be used to minimize functions, and the solution type has a continuous interpolation which can be used. Let's look for the local optima for the 4th variable around `t=20`. Thus, our optimization function is:

```@example minmax
f(t, _) = sol(first(t), idxs = 4)
```

`first(t)` is the same as `t[1]` which transforms the array of size 1 into a number. `idxs=4` is the same as `sol(first(t))[4]` but does the calculation without a temporary array and thus is faster. To find a local minimum, we can solve the optimization problem where the loss
function is `f`:

```@example minmax
import Optimization as OPT, OptimizationNLopt as OptNL, ForwardDiff
optf = OPT.OptimizationFunction(f, OPT.AutoForwardDiff())
min_guess = 18.0
optprob = OPT.OptimizationProblem(optf, [min_guess], lb = [0.0], ub = [100.0])
opt = OPT.solve(optprob, OptNL.NLopt.LD_LBFGS())
```

From this printout, we see that the minimum is at `t=18.63` and the value is `-2.79e-2`. We
can get these in code-form via:

```@example minmax
println(opt.u)
```

To get the maximum, we just minimize the negative of the function:

```@example minmax
fminus(t, _) = -sol(first(t), idxs = 4)

optf = OPT.OptimizationFunction(fminus, OPT.AutoForwardDiff())
min_guess = 22.0
optprob2 = OPT.OptimizationProblem(optf, [min_guess], lb = [0.0], ub = [100.0])
opt2 = OPT.solve(optprob2, OptNL.NLopt.LD_LBFGS())
```

Let's add the maxima and minima to the plots:

```@example minmax
Plots.plot(sol, vars = (0, 4), plotdensity = 10000)
Plots.scatter!([opt.u], [opt.minimum], label = "Local Min")
Plots.scatter!([opt2.u], [-opt2.minimum], label = "Local Max")
```

### Global Optimization

If we instead want to find global maxima and minima, we need to look somewhere else.
There are many choices for this. A pure Julia option are the
[BlackBoxOptim solvers within Optimization.jl](https://docs.sciml.ai/Optimization/stable/optimization_packages/blackboxoptim/),
but I will continue the story with the  OptimizationNLopt methods. To do this, we simply
swap out to one of the
[global optimizers in the list](https://docs.sciml.ai/Optimization/stable/optimization_packages/nlopt/).
Let's try `GN_ORIG_DIRECT_L`:

```@example minmax
gopt = OPT.solve(optprob, OptNL.NLopt.GN_ORIG_DIRECT_L())
gopt2 = OPT.solve(optprob2, OptNL.NLopt.GN_ORIG_DIRECT_L())

@show gopt.u, gopt2.u
```

```@example minmax
Plots.plot(sol, vars = (0, 4), plotdensity = 10000)
Plots.scatter!([gopt.u], [gopt.minimum], label = "Global Min")
Plots.scatter!([gopt2.u], [-gopt2.minimum], label = "Global Max")
```
