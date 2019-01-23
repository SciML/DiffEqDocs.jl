# Steady State Solvers

`solve(prob::SteadyStateProblem,alg;kwargs)`

Solves for the steady states in the problem defined by `prob` using the algorithm
`alg`. If no algorithm is given, a default algorithm will be chosen.

## Recommended Methods

`DynamicSS` is a good choice if you think you may have multiple steady states
or a bad initial guess. `SSRootfind` can be faster if you have a good initial
guess. For `DynamicSS`, in many cases an adaptive stiff solver, like a
Rosenbrock method (`Rodas5` or `CVODE_BDF`), is a good way to allow for very
large time steps as the steady state approaches. Note that if you use `CVODE_BDF`
you may need to give a starting `dt` via `dt=....`.

## Full List of Methods

### SteadyStateDiffEq.jl

- `SSRootfind` : Uses a rootfinding algorithm to find a steady state. Defaults
  to using NLsolve.jl. A different algorithm can be specified via the `nlsolve`
  keyword argument.
- `DynamicSS` : Uses an ODE solver to find the steady state. Automatically
  terminates when close to the steady state.
  `DynamicSS(alg;abstol=1e-8,reltol=1e-6,tspan=Inf)` requires that an
  ODE algorithm is given as the first argument.  The absolute and
  relative tolerances specify the termination conditions on the
  derivative's closeness to zero.  This internally uses the
  `TerminateSteadyState` callback from the Callback Library.  The
  simulated time for which given ODE is solved can be limited by
  `tspan`.  If `tspan` is a number, it is equivalent to passing
  `(zero(tspan), tspan)`.

Example usage:

```julia
sol = solve(prob,SSRootfind())
sol = solve(prob,DynamicSS(Tsit5()))
sol = solve(prob,DynamicSS(CVODE_BDF()),dt=1.0)
```
