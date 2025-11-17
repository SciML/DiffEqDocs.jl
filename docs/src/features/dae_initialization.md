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
DiffEqBase.BrownFullBasicInit
DiffEqBase.ShampineCollocationInit
```

## ⚠️ WARNING: NoInit Usage

!!! warn "Use NoInit at your own risk"
    **`NoInit()` should almost never be used.** No algorithm has any guarantee of correctness if `NoInit()` is used with inconsistent initial conditions. Users should almost always use `CheckInit()` instead for safety.

    **Important:**
    - Any issues opened that are using `NoInit()` will be immediately closed
    - Allowing incorrect initializations is not a supported part of the system
    - Using `NoInit()` with inconsistent conditions can lead to:
      + Solver instability and crashes
      + Incorrect results that may appear plausible
      + Undefined behavior in the numerical algorithms
      + Silent corruption of the solution

    **When to use `CheckInit()` instead:**
    - When you believe your initial conditions are consistent
    - When you want to skip automatic modification of initial conditions
    - When you need to verify your manual initialization

    The only valid use case for `NoInit()` is when you are 100% certain your conditions are consistent AND you need to skip the computational cost of verification for performance reasons in production code that has been thoroughly tested.

## Algorithm Selection Guide

| Algorithm                   | When to Use                                               | Modifies Variables        |
|:--------------------------- |:--------------------------------------------------------- |:------------------------- |
| `DefaultInit()`             | Default choice - automatically selects appropriate method | Depends on selection      |
| `CheckInit()`               | When you've computed consistent conditions yourself       | No (verification only)    |
| `NoInit()`                  | ⚠️ **AVOID** - Only for verified consistent conditions    | No                        |
| `OverrideInit()`            | With ModelingToolkit problems                             | Yes (uses custom problem) |
| `BrownFullBasicInit()`      | For index-1 DAEs with `differential_vars`                 | Algebraic variables only  |
| `ShampineCollocationInit()` | For general DAEs without structure information            | All variables             |

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

# BrownFullBasicInit will fix the inconsistent du0
sol = solve(prob, DFBDF(), initializealg = BrownFullBasicInit())
```

### Example 2: Checking Consistency (Recommended over NoInit)

```julia
# If you've computed consistent conditions yourself
u0_consistent = [1.0, 0.0, 0.0]
du0_consistent = [0.0, -1.0, compute_tension(u0_consistent, p)]

prob2 = DAEProblem(pendulum!, du0_consistent, u0_consistent, tspan, p,
    differential_vars = [true, true, false])

# RECOMMENDED: Verify they're consistent with CheckInit
sol = solve(prob2, DFBDF(), initializealg = CheckInit())

# NOT RECOMMENDED: NoInit skips all checks - use at your own risk!
# sol = solve(prob2, DFBDF(), initializealg = NoInit())  # ⚠️ DANGEROUS
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
sys = structural_simplify(pendulum)

# ModelingToolkit provides initialization_data
prob = DAEProblem(sys, [x => 1.0, y => 0.0], (0.0, 10.0), [g => 9.81, L => 1.0])

# DefaultInit will use OverrideInit with ModelingToolkit's initialization_data
sol = solve(prob, DFBDF())  # Automatic initialization!
```

## Algorithm-Specific Options

Both OrdinaryDiffEq and Sundials support the same initialization algorithms through the `initializealg` parameter:

### OrdinaryDiffEq and Sundials

```julia
using OrdinaryDiffEq
# or
using Sundials

# Use Brown's algorithm to fix inconsistent conditions
sol = solve(prob, DFBDF(), initializealg = BrownFullBasicInit())  # OrdinaryDiffEq
sol = solve(prob, IDA(), initializealg = BrownFullBasicInit())    # Sundials

# Use Shampine's collocation method for general DAEs
sol = solve(prob, DFBDF(), initializealg = ShampineCollocationInit())  # OrdinaryDiffEq
sol = solve(prob, IDA(), initializealg = ShampineCollocationInit())    # Sundials

# RECOMMENDED: Verify conditions are consistent
sol = solve(prob, DFBDF(), initializealg = CheckInit())  # OrdinaryDiffEq
sol = solve(prob, IDA(), initializealg = CheckInit())    # Sundials

# NOT RECOMMENDED: Skip all initialization checks (dangerous!)
# sol = solve(prob, DFBDF(), initializealg = NoInit())  # ⚠️ USE AT YOUR OWN RISK
# sol = solve(prob, IDA(), initializealg = NoInit())    # ⚠️ USE AT YOUR OWN RISK
```

## Troubleshooting

### Common Issues and Solutions

1. **"Initial conditions are not consistent" error**
   + Ensure your `du0` satisfies the DAE constraints at `t0`
   + Try using `BrownFullBasicInit()` or `ShampineCollocationInit()` instead of `CheckInit()`
   + Check that `differential_vars` correctly identifies differential vs algebraic variables

2. **Initialization fails to converge**
   + Relax tolerances if using extended versions
   + Try a different initialization algorithm
   + Provide a better initial guess for algebraic variables
   + **Check if your DAE is index-1**: The system may be higher-index (see below)

3. **Solver fails immediately after initialization**
   + The initialization might have found a consistent but numerically unstable point
   + Try tightening initialization tolerances
   + Check problem scaling and consider non-dimensionalization

4. **DAE is not index-1 (higher-index DAE)**
   + Many initialization algorithms only work reliably for index-1 DAEs
   + **To check if your DAE is index-1**: The Jacobian of the algebraic equations with respect to the algebraic variables must be non-singular
   + **Solution**: Use ModelingToolkit.jl to analyze and potentially reduce the index:
     ```julia
     using ModelingToolkit

     # Define your system with ModelingToolkit
     @named sys = ODESystem(eqs, t, vars, params)

     # Analyze and reduce the index (structural_simplify handles this in v10+)
     sys_reduced = structural_simplify(sys)

     # The reduced system will be index-1 and easier to initialize
     prob = DAEProblem(sys_reduced, [], (0.0, 10.0), params)
     ```
   + ModelingToolkit can automatically detect the index and apply appropriate transformations
   + After index reduction, standard initialization algorithms will work more reliably

## Performance Tips

 1. **Use `differential_vars` when possible**: This helps initialization algorithms understand problem structure
 2. **Provide good initial guesses**: Even when using automatic initialization, starting closer to the solution helps
 3. **Consider problem-specific initialization**: For complex systems, custom initialization procedures may be more efficient
 4. **Use `CheckInit()` when appropriate**: If you know conditions are consistent, skip unnecessary computation

## References

  - Brown, P. N., Hindmarsh, A. C., & Petzold, L. R. (1998). Consistent initial condition calculation for differential-algebraic systems. SIAM Journal on Scientific Computing, 19(5), 1495-1512.
  - Shampine, L. F. (2002). Consistent initial condition for differential-algebraic systems. SIAM Journal on Scientific Computing, 22(6), 2007-2026.
