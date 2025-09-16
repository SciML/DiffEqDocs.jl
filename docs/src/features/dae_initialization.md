# DAE Initialization

DAE (Differential-Algebraic Equation) problems often require special initialization procedures to ensure that the initial conditions are consistent with the algebraic constraints. The DifferentialEquations.jl ecosystem provides several initialization algorithms to handle this automatically or to verify that your provided initial conditions are already consistent.

## The Initialization Problem

DAEs have the general form:

```math
M \frac{du}{dt} = f(u, p, t)
```

where `M` is a (possibly singular) mass matrix. For the initial conditions `u₀` and `du₀` to be consistent, they must satisfy:

```math
f(du₀, u₀, p, t₀) = 0
```

for fully implicit DAEs, or the equivalent constraint for semi-explicit DAEs. Finding consistent initial conditions is a nonlinear problem that must be solved before time integration can begin.

## Available Initialization Algorithms

The `initializealg` keyword argument to `solve` controls how initialization is performed. All algorithms are documented with their docstrings:

```@docs
DiffEqBase.DefaultInit
DiffEqBase.CheckInit
DiffEqBase.NoInit
DiffEqBase.OverrideInit
DiffEqBase.BrownBasicInit
DiffEqBase.ShampineCollocationInit
```

## Algorithm Selection Guide

| Algorithm | When to Use | Modifies Variables |
|-----------|-------------|-------------------|
| `DefaultInit()` | Default choice - automatically selects appropriate method | Depends on selection |
| `CheckInit()` | When you've computed consistent conditions yourself | No (verification only) |
| `NoInit()` | When conditions are known to be perfectly consistent | No |
| `OverrideInit()` | With ModelingToolkit problems | Yes (uses custom problem) |
| `BrownBasicInit()` | For index-1 DAEs with `differential_vars` | Algebraic variables only |
| `ShampineCollocationInit()` | For general DAEs without structure information | All variables |

## Examples

### Example 1: Simple Pendulum DAE

```julia
using DifferentialEquations

function pendulum!(res, du, u, p, t)
    x, y, T = u
    dx, dy, dT = du
    g, L = p

    res[1] = dx - du[1]
    res[2] = dy - du[2]
    res[3] = x^2 + y^2 - L^2  # Algebraic constraint
end

u0 = [1.0, 0.0, 0.0]  # Initial position
du0 = [0.0, 0.0, 0.0]  # Initial velocity (inconsistent!)
p = [9.81, 1.0]  # g, L
tspan = (0.0, 10.0)

prob = DAEProblem(pendulum!, du0, u0, tspan, p,
                  differential_vars = [true, true, false])

# BrownBasicInit will fix the inconsistent du0
sol = solve(prob, DFBDF(), initializealg = BrownBasicInit())
```

### Example 2: Checking Consistency

```julia
# If you've computed consistent conditions yourself
u0_consistent = [1.0, 0.0, 0.0]
du0_consistent = [0.0, -1.0, compute_tension(u0_consistent, p)]

prob2 = DAEProblem(pendulum!, du0_consistent, u0_consistent, tspan, p,
                   differential_vars = [true, true, false])

# Just verify they're consistent
sol = solve(prob2, DFBDF(), initializealg = CheckInit())
```

### Example 3: ModelingToolkit Integration

When using ModelingToolkit, initialization information is often included automatically:

```julia
using ModelingToolkit, DifferentialEquations

@variables t x(t) y(t) T(t)
@parameters g L
D = Differential(t)

eqs = [
    D(x) ~ -T * x/L,
    D(y) ~ -T * y/L - g,
    x^2 + y^2 ~ L^2
]

@named pendulum = ODESystem(eqs, t, [x, y, T], [g, L])
sys = structural_simplify(dae_index_lowering(pendulum))

# ModelingToolkit provides initialization_data
prob = DAEProblem(sys, [x => 1.0, y => 0.0], (0.0, 10.0), [g => 9.81, L => 1.0])

# DefaultInit will use OverrideInit with ModelingToolkit's initialization_data
sol = solve(prob, DFBDF())  # Automatic initialization!
```

## Algorithm-Specific Options

Some solvers provide extended versions of these algorithms with additional options:

### OrdinaryDiffEq Extensions

OrdinaryDiffEq provides extended versions with additional parameters:

```julia
using OrdinaryDiffEq

# Shampine with custom initial dt and nonlinear solver
sol = solve(prob, DFBDF(),
            initializealg = OrdinaryDiffEqCore.ShampineCollocationInitExt(initdt = 0.001))

# Brown with custom absolute tolerance
sol = solve(prob, DFBDF(),
            initializealg = OrdinaryDiffEqCore.BrownFullBasicInit(abstol = 1e-10))
```

### Sundials (IDA)

The IDA solver from Sundials.jl uses initialization through the `initializealg` parameter:

```julia
using Sundials

# Use Brown's algorithm
sol = solve(prob, IDA(), initializealg = BrownBasicInit())

# Skip initialization if you know conditions are consistent
sol = solve(prob, IDA(), initializealg = NoInit())
```

## Troubleshooting

### Common Issues and Solutions

1. **"Initial conditions are not consistent" error**
   - Ensure your `du0` satisfies the DAE constraints at `t0`
   - Try using `BrownBasicInit()` or `ShampineCollocationInit()` instead of `CheckInit()`
   - Check that `differential_vars` correctly identifies differential vs algebraic variables

2. **Initialization fails to converge**
   - Relax tolerances if using extended versions
   - Try a different initialization algorithm
   - Provide a better initial guess for algebraic variables

3. **Solver fails immediately after initialization**
   - The initialization might have found a consistent but numerically unstable point
   - Try tightening initialization tolerances
   - Check problem scaling and consider non-dimensionalization

## Performance Tips

1. **Use `differential_vars` when possible**: This helps initialization algorithms understand problem structure
2. **Provide good initial guesses**: Even when using automatic initialization, starting closer to the solution helps
3. **Consider problem-specific initialization**: For complex systems, custom initialization procedures may be more efficient
4. **Use `CheckInit()` when appropriate**: If you know conditions are consistent, skip unnecessary computation

## References

- Brown, P. N., Hindmarsh, A. C., & Petzold, L. R. (1998). Consistent initial condition calculation for differential-algebraic systems. SIAM Journal on Scientific Computing, 19(5), 1495-1512.
- Shampine, L. F. (2002). Consistent initial condition for differential-algebraic systems. SIAM Journal on Scientific Computing, 22(6), 2007-2026.