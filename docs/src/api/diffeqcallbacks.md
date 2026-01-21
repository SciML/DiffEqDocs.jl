# [DiffEqCallbacks.jl](@id diffeqcallbacks_api)

DiffEqCallbacks.jl provides a library of pre-built callbacks for use with the
SciML differential equation solvers. These include saving callbacks, manifold
projection, domain constraints, and more.

## Installation

DiffEqCallbacks.jl is included with DifferentialEquations.jl. To use it standalone:

```julia
using Pkg
Pkg.add("DiffEqCallbacks")
import DiffEqCallbacks
```

## Callback APIs

### Manifold Projection Callbacks

```@docs
DiffEqCallbacks.ManifoldProjection
```

### Saving Callbacks

```@docs
DiffEqCallbacks.SavingCallback
DiffEqCallbacks.SavedValues
```

### Domain Callbacks

```@docs
DiffEqCallbacks.PositiveDomain
DiffEqCallbacks.GeneralDomain
```

### Stepping Callbacks

```@docs
DiffEqCallbacks.StepsizeLimiter
DiffEqCallbacks.FunctionCallingCallback
```

### Termination Callbacks

```@docs
DiffEqCallbacks.TerminateSteadyState
```

### Iterative Callbacks

```@docs
DiffEqCallbacks.IterativeCallback
DiffEqCallbacks.PeriodicCallback
```

### Preset Time Callbacks

```@docs
DiffEqCallbacks.PresetTimeCallback
```

### AutoAbstol

```@docs
DiffEqCallbacks.AutoAbstol
```
