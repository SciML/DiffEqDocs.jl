# [SimpleDiffEq.jl](@id simplediffeq_api)

SimpleDiffEq.jl provides simplified implementations of a few ODE and SDE solvers.
They are primarily designed for experimentation and offer shorter compile times.
They have limitations compared to OrdinaryDiffEq.jl and StochasticDiffEq.jl and
are generally not faster.

Note that this package is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use SimpleDiffEq.jl:

```julia
using Pkg
Pkg.add("SimpleDiffEq")
import SimpleDiffEq
```

## ODE Solvers

The following ODE solvers are available:

  - ``SimpleTsit5`` - A fixed timestep integrator form of Tsit5. Not compatible with events.
  - ``SimpleATsit5`` - An adaptive Tsit5 with an interpolation in its simplest form. Not compatible with events.
  - ``GPUSimpleATsit5`` - A version of SimpleATsit5 without the integrator interface. Only allows `solve`.
  - ``SimpleEuler`` - A fixed timestep bare-bones Euler implementation with integrators.
  - ``SimpleRK4`` - A fixed timestep bare-bones RK4 implementation with integrators.
  - ``GPUSimpleVern7`` - A fully static Vern7 for specialized compilation to accelerators like GPUs and TPUs.
  - ``GPUSimpleVern9`` - A fully static Vern9 for specialized compilation to accelerators like GPUs and TPUs.
  - ``GPUSimpleTsit5`` - A fully static Tsit5 for specialized compilation to accelerators.
  - ``GPUSimpleRK4`` - A fully static RK4 for specialized compilation to accelerators.
  - ``GPUSimpleEuler`` - A fully static Euler for specialized compilation to accelerators.
  - ``LoopEuler`` - A fixed timestep bare-bones Euler. Not compatible with events or the integrator interface.
  - ``LoopRK4`` - A fixed timestep bare-bones RK4. Not compatible with events or the integrator interface.

## SDE Solvers

The following SDE solvers are available:

  - ``SimpleEM`` - A fixed timestep solve method for Euler-Maruyama. Only works with non-colored Gaussian noise.
