# Verbosity Specification with SciMLLogging.jl

DifferentialEquations.jl uses [SciMLLogging.jl](https://docs.sciml.ai/SciMLLogging/stable/) to provide users with fine-grained control over logging and diagnostic output during ODE solving. The `DEVerbosity` struct allows you to customize which messages are displayed, from critical errors to detailed debugging information.

## Basic Usage

Pass a `DEVerbosity` object to `solve` or `init` using the `verbose` keyword argument:

```julia
using OrdinaryDiffEq

# Define an ODE problem
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8/3) * u[3]
end

u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)

# Solve with detailed verbosity
verbose = DEVerbosity(SciMLLogging.Detailed())
sol = solve(prob, Tsit5(), verbose = verbose)
```

## Example Use Cases

### Debugging Algorithm Switching

When using auto-switching algorithms, see when and why switches occur:

```julia
using ODEProblemLibrary: prob_ode_vanderpol_stiff

verbose = DEVerbosity(alg_switch = SciMLLogging.InfoLevel())
sol = solve(prob_ode_vanderpol_stiff, AutoTsit5(Rodas5()), verbose = verbose)
```

### Monitoring Stiff Solver Performance

Get detailed information about Jacobian updates, factorizations, and linear solver behavior:

```julia
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
end

u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 0.1)
p = [0.04, 3e7, 1e4]
prob = ODEProblem(rober, u0, tspan, p)

# Monitor performance-related events and linear solver details
verbose = DEVerbosity(
    performance = SciMLLogging.InfoLevel(),
    linear_verbosity = SciMLLogging.Detailed()
)
sol = solve(prob, Rosenbrock23(), verbose = verbose)
```

### Debugging Step Acceptance/Rejection

For adaptive methods, monitor which steps are accepted or rejected:

```julia
verbose = DEVerbosity(
    step_accepted = SciMLLogging.InfoLevel(),
    step_rejected = SciMLLogging.InfoLevel()
)
sol = solve(prob, Tsit5(), verbose = verbose)
```

### Silent Operation

Completely disable all output:

```julia
verbose = DEVerbosity(SciMLLogging.None())
solve(prob, Tsit5(), verbose = verbose)
```

## Using with `init`

The verbosity settings work the same way with `init` for manual stepping:

```julia
verbose = DEVerbosity(
    alg_switch = SciMLLogging.InfoLevel(),
    linear_verbosity = SciMLLogging.Detailed()
)
integrator = init(prob, Rosenbrock23(), verbose = verbose, dt = 1e-3)

# Step through the solution
step!(integrator)
```

## Controlling Subsolver Verbosity

When using implicit methods, `linear_verbosity` and `nonlinear_verbosity` are automatically passed to the linear and nonlinear solver caches, enabling diagnostics from LinearSolve.jl and NonlinearSolve.jl.

### Using SciMLLogging Presets

The simplest approach is to pass a SciMLLogging preset, which is automatically converted to the appropriate verbosity type:

```julia
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg

# Enable detailed nonlinear solver diagnostics
verbose = DEVerbosity(nonlinear_verbosity = SciMLLogging.Detailed())
sol = solve(prob, ImplicitEuler(nlsolve = NonlinearSolveAlg()), verbose = verbose)
```

### Using LinearVerbosity and NonlinearVerbosity Explicitly

You can also explicitly use `LinearVerbosity` and `NonlinearVerbosity` types from LinearSolve.jl and NonlinearSolve.jl:

```julia
using LinearSolve: LinearVerbosity
using NonlinearSolve: NonlinearVerbosity

# Create explicit verbosity objects
linear_verb = LinearVerbosity(SciMLLogging.Detailed())
nonlinear_verb = NonlinearVerbosity(SciMLLogging.Minimal())

# Pass them to DEVerbosity
verbose = DEVerbosity(
    linear_verbosity = linear_verb,
    nonlinear_verbosity = nonlinear_verb
)

# Use with a stiff solver
sol = solve(prob, Rosenbrock23(), verbose = verbose)
```

## Combining Settings Flexibly

Individual field settings override group settings, allowing precise control:

```julia
# Set all error_control messages to WarnLevel,
# but make step_rejected silent and dt_NaN an error
verbose = DEVerbosity(
    error_control = SciMLLogging.WarnLevel(),
    step_rejected = SciMLLogging.Silent(),
    dt_NaN = SciMLLogging.ErrorLevel()
)
sol = solve(prob, Tsit5(), verbose = verbose)
```

## API Reference

```@docs
DEVerbosity
```