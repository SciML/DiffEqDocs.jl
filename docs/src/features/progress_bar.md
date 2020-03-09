# Progress Bar Integration

DifferentialEquations.jl integrates with the progress bars in order to make
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
the timesteps may (and will) change. In Juno, scrolling over the progressbar one will
also see the current timestep. This can be used to track the solution's progress
and find tough locations for the solvers.
