# Callback Library

DiffEqCallbackLibrary.jl provides a library of various helpful callbacks which
can be used with any component solver which implements the callback interface.
As the package is currently unregistered, it must first be installed via:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqCallbacks.jl")
using DiffEqCallbacks
```

## Callbacks

### AutoAbstol

Many problem solving environments [such as MATLAB](https://www.mathworks.com/help/simulink/gui/absolute-tolerance.html)
provide a way to automatically adapt the absolute tolerance to the problem. This
helps the solvers automatically "learn" what appropriate limits are. Via the
callback interface, DiffEqCallbacks.jl implements a callback `AutoAbstol` which
has the same behavior as the MATLAB implementation, that is the absolute tolerance
starts at `init_curmax` (default `1-e6`), and at each iteration it is set
to the maximum value that the state has thus far reached times the relative tolerance.

To generate the callback, use the constructor:

```julia
AutoAbstol(save=true;init_curmax=1e-6)
```
