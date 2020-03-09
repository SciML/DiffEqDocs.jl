# Progress Bar Integration

DifferentialEquations.jl integrates with the Juno progress bar in order to make
long calculations more manageable. By default this feature is off for ODE and
SDE solvers, but can be turned on via the keyword argument `progressbar=true`.
The progress bar updates every `progress_steps` timesteps, which has a default
value of 1000. Note that making this value really low could cause a performance
hit, though from some basic testing it seems that with updates of at least
1000 steps on number (the fastest problems) there's no discernable performance degradation,
giving a high upper bound.

Note that the progressbar also includes a time estimate. This time-estimate is provided
by linear extrapolation for how long it has taken to get to what percentage. For
adaptive timestepping methods this should only be used as a rough estimate since
the timesteps may (and will) change. By scrolling over the progressbar one will
also see the current timestep. This can be used to track the solution's progress
and find tough locations for the solvers.

## Using Progress Bars Outside Juno

To use the progress bars outside of Juno, use [TerminalLoggers.jl](https://github.com/c42f/TerminalLoggers.jl).
The following is an example for redirecting the logging to the terminal:

```julia
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using OrdinaryDiffEq

solve(
    ODEProblem((u, p, t) -> (sleep(0.01); -u), 1.0, nothing),
    Euler();
    dt = 0.5,
    tspan = (0.0, 1000.0),
    progress = true,
    progress_steps = 1,
)
```

To do this by default, [follow these direction to add TerminalLogging to your startup.jl](https://c42f.github.io/TerminalLoggers.jl/stable/#Installation-and-setup-1).
