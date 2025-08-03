# Progress Bar Integration

DifferentialEquations.jl integrates with the Juno progress bar in order to make
long calculations more manageable. By default, this feature is off for ODE and
SDE solvers, but can be turned on via the keyword argument `progress=true`.
The progress bar updates every `progress_steps` timesteps, which has a default
value of 1000. Note that making this value really low could cause a performance
hit, though from some basic testing it seems that with updates of at least
1000 steps on number (the fastest problems), there's no discernible performance degradation,
giving a high upper bound.

Note that the progress bar also includes a time estimate. This time-estimate is provided
by linear extrapolation for how long it has taken to get to what percentage. For
adaptive timestepping methods this should only be used as a rough estimate since
the timesteps may (and will) change. By scrolling over the progress bar, one will
also see the current timestep. This can be used to track the solution's progress
and find tough locations for the solvers.

## Using Progress Bars with VS Code

If using VS Code, progress bars are enabled via the [ProgressLogging.jl](https://github.com/JuliaLogging/ProgressLogging.jl) package.
For example:

```julia
import OrdinaryDiffEq as ODE, ProgressLogging
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1000000.0)
prob = ODE.ODEProblem(lorenz!, u0, tspan)
sol = ODE.solve(prob, ODE.Tsit5(), progress = true)
```

## Using Progress Bars in the Terminal

```julia
import OrdinaryDiffEq as ODE, TerminalLoggers
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1000000.0)
prob = ODE.ODEProblem(lorenz!, u0, tspan)
sol = ODE.solve(prob, ODE.Tsit5(), progress = true)
```

To use progress bars in the terminal, use [TerminalLoggers.jl](https://github.com/JuliaLogging/TerminalLoggers.jl).
[Follow these directions to add TerminalLogging to your startup.jl](https://julialogging.github.io/TerminalLoggers.jl/stable/#Installation-and-setup-1),
if you want it enabled by default.

Otherwise, follow the example down below. Note that `global_logger` is initialized
before any other Julia call. This step is crucial. Otherwise, no logging will
appear in the terminal.

```julia
import OrdinaryDiffEq as ODE
import Logging: global_logger
import TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1000000.0)
prob = ODE.ODEProblem(lorenz!, u0, tspan)
sol = ODE.solve(prob, ODE.Tsit5(), progress = true, maxiters = 1e8)
```
