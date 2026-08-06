# [DASSL.jl](@id dassl_api)

DASSL.jl is a native Julia implementation of the DASSL algorithm for solving
differential-algebraic equations (DAEs) within the SciML interface.

Note that this package is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use DASSL.jl:

```julia
using Pkg
Pkg.add("DASSL")
import DASSL
```

These methods can be used independently of the rest of DifferentialEquations.jl.

## DAE Solvers

The following DAE solver is available:

  - `dassl` - A native Julia implementation of the DASSL algorithm. This is an
    adaptive-order, adaptive-stepsize backward differentiation formula (BDF)
    method for solving stiff systems of differential-algebraic equations.

### Example Usage

```julia
using DASSL, DiffEqBase

# Define the residual function for the DAE
function f!(res, du, u, p, t)
    res[1] = du[1] - u[2]
    res[2] = u[1] + u[2] - 1.0
end

# Initial conditions
u0 = [1.0, 0.0]
du0 = [0.0, 1.0]
tspan = (0.0, 10.0)

# Create and solve the problem
prob = DAEProblem(f!, du0, u0, tspan)
sol = solve(prob, dassl())
```
