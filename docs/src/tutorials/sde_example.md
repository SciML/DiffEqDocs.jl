# Stochastic Differential Equations

!!! note
    
    This tutorial assumes you have read the [Ordinary Differential Equations tutorial](@ref ode_example).

## Example 1: Scalar SDEs

In this example, we will solve the equation

```math
du = f(u,p,t)dt + g(u,p,t)dW,
```

where ``f(u,p,t)=αu`` and ``g(u,p,t)=βu``. We know via Stochastic calculus that the
solution to this equation is

```math
u(t,Wₜ)=u₀\exp\left[\left(α-\frac{β^2}{2}\right)t+βWₜ\right].
```

To solve this numerically, we define a stochastic problem type using `SDEProblem` by specifying `f(u, p, t)`, `g(u, p, t)`, and the initial condition:

```@example sde
import DifferentialEquations as DE
α = 1
β = 1
u₀ = 1 / 2
f(u, p, t) = α * u
g(u, p, t) = β * u
dt = 1 // 2^(4)
tspan = (0.0, 1.0)
prob = DE.SDEProblem(f, g, u₀, tspan)
```

The `solve` interface is then the same as ODEs. Here, we will use the classic
Euler-Maruyama algorithm `EM` and plot the solution:

```@example sde
sol = DE.solve(prob, DE.EM(), dt = dt)
import Plots
Plots.plot(sol)
```

### Using Higher Order Methods

One unique feature of DifferentialEquations.jl is that higher-order methods for
stochastic differential equations are included. To illustrate it, let us compare the
accuracy of the `DE.EM()` method and a higher-order method `DE.SRIW1()` with the analytical solution.
This is a good way to judge the accuracy of a given algorithm, and is also useful
to test convergence of new methods being developed. To setup our problem, we define
`u_analytic(u₀, p, t, W)` and pass it to the `SDEFunction` as:

```@example sde
u_analytic(u₀, p, t, W) = u₀ * exp((α - (β^2) / 2) * t + β * W)
ff = DE.SDEFunction(f, g, analytic = u_analytic)
prob = DE.SDEProblem(ff, u₀, (0.0, 1.0))
```

We can now compare the `DE.EM()` solution with the analytic one:

```@example sde
sol = DE.solve(prob, DE.EM(), dt = dt)
Plots.plot(sol, plot_analytic = true)
```

Now, we choose a higher-order solver `DE.SRIW1()` for better accuracy. By default,
the higher order methods are adaptive. Let's first switch off adaptivity and
compare the numerical and analytic solutions :

```@example sde
sol = DE.solve(prob, DE.SRIW1(), dt = dt, adaptive = false)
Plots.plot(sol, plot_analytic = true)
```

Now, let's allow the solver to automatically determine a starting `dt`. This estimate
at the beginning is conservative (small) to ensure accuracy.

```@example sde
sol = DE.solve(prob, DE.SRIW1())
Plots.plot(sol, plot_analytic = true)
```

We can instead start the method with a larger `dt` by passing it to `solve`:

```@example sde
sol = DE.solve(prob, DE.SRIW1(), dt = dt)
Plots.plot(sol, plot_analytic = true)
```

### Ensemble Simulations

Instead of solving single trajectories, we can turn our problem into a `EnsembleProblem`
to solve many trajectories all at once. This is done by the `EnsembleProblem`
constructor:

```@example sde
ensembleprob = DE.EnsembleProblem(prob)
```

The solver commands are defined [at the Parallel Ensemble Simulations page](@ref ensemble).
For example, we can choose to have 1000 trajectories via `trajectories=1000`. In addition,
this will automatically parallelize using Julia native parallelism if extra processes
are added via `addprocs()`, but we can change this to use multithreading via
`DE.EnsembleThreads()`. Together, this looks like:

```@example sde
sol = DE.solve(ensembleprob, DE.EnsembleThreads(), trajectories = 1000)
```

!!! warn
    
    If you use a custom noise process, you might need to specify it in a custom prob_func
    in the EnsembleProblem constructor, as each trajectory needs its own noise process.

Many more controls are defined at the [Ensemble simulations page](@ref ensemble),
including analysis tools.
A very simple analysis can be done with the `EnsembleSummary`, which builds
mean/var statistics and has an associated plot recipe. For example, we can get
the statistics at every `0.01` timesteps and plot the average + error using:

```@example sde
import DifferentialEquations as DE
summ = DE.EnsembleSummary(sol, 0:0.01:1)
Plots.plot(summ, labels = "Middle 95%")
summ = DE.EnsembleSummary(sol, 0:0.01:1; quantiles = [0.25, 0.75])
Plots.plot!(summ, labels = "Middle 50%", legend = true)
```

Additionally, we can easily calculate the correlation between the values at `t=0.2`
and `t=0.7` via

```@example sde
DE.timepoint_meancor(sol, 0.2, 0.7) # Gives both means and then the correlation coefficient
```

## Example 2: Systems of SDEs with Diagonal Noise

In general, a system of SDEs

```math
du = f(u,p,t)dt + g(u,p,t)dW,
```

where `u` is now a vector of variables, `f` is a vector, and `g` is a matrix, is numerically integrated in the same way as ODEs. A common scenario, which is the default for DifferentialEquations.jl, is that every variable in the system gets a different random kick. This is the case when `g` is a diagonal matrix. Correspondingly, we say that we have a diagonal noise. 

We handle this in a simple manner by defining the deterministic part `f!(du,u,p,t)` and the stochastic part
`g!(du2,u,p,t)` as in-place functions, but note that our convention is that the function `g!` only defines and modifies the diagonal entries of `g` matrix.

As an example, we consider a stochastic variant of the Lorenz equations, where we add to each of `u₁, u₂, u₃` their own simple noise `3*N(0,dt)`. Here, `N` is the normal distribution and `dt` is the time step. This is done as follows:

```@example sde2
import DifferentialEquations as DE
import Plots

function f!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

function g!(du, u, p, t)  # It actually represents a diagonal matrix [3.0 0 0; 0 3.0 0; 0 0 3.0]
    du[1] = 3.0
    du[2] = 3.0
    du[3] = 3.0
end

prob_sde_lorenz = SDEProblem(f!, g!, [1.0, 0.0, 0.0], (0.0, 10.0))
sol = DE.solve(prob_sde_lorenz)
Plots.plot(sol, idxs = (1, 2, 3))
```

Note that it's okay for the noise function to mix terms. For example

```@example sde2
function g!(du, u, p, t)
    du[1] = sin(u[3]) * 3.0
    du[2] = u[2] * u[1] * 3.0
    du[3] = 3.0
end
```

is a valid noise function.

## Example 3: Systems of SDEs with Scalar Noise

In this example, we'll solve a system of SDEs with scalar noise. This means that
the same noise process is applied to all SDEs. We need to define a
scalar noise process
[using the Noise Process interface](@ref noise_process).
Since we want a `WienerProcess` that starts at `0.0` at time `0.0`, we use the
command `W = DE.WienerProcess(0.0,0.0,0.0)` to define the Brownian motion we want,
and then give this to the `noise` option in the `SDEProblem`. For a full example,
let's solve a linear SDE with scalar noise using a high order algorithm:

```@example sde3
import DifferentialEquations as DE
import Plots
f!(du, u, p, t) = (du .= u)
g!(du, u, p, t) = (du .= u)
u0 = rand(4, 2)

W = DE.WienerProcess(0.0, 0.0, 0.0)
prob = DE.SDEProblem(f!, g!, u0, (0.0, 1.0), noise = W)
sol = DE.solve(prob, DE.SRIW1())
Plots.plot(sol)
```

## Example 4: Systems of SDEs with Non-Diagonal Noise

In the previous examples we had diagonal noise, that is a vector of random numbers
`dW` whose size matches the output of `g` where the noise is applied element-wise,
and scalar noise where a single random variable is applied to all dependent variables.
However, a more general type of noise allows for the terms to linearly mixed via `g`
being a general nondiagonal matrix.

Note that nonlinear mixings are not SDEs but fall under the more general class of
random ordinary differential equations (RODEs) which have a
[separate set of solvers](@ref rode_example).

Let's define a problem with four Wiener processes and two dependent random variables.
In this case, we will want the output of `g` to be a 2x4 matrix, such that the solution
is `g(u,p,t)*dW`, the matrix multiplication. For example, we can do the following:

```@example sde4
import DifferentialEquations as DE
f!(du, u, p, t) = du .= 1.01u
function g!(du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[1]
    du[2, 1] = 1.2u[2]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    du[2, 4] = 1.8u[2]
end
prob = DE.SDEProblem(f!, g!, ones(2), (0.0, 1.0), noise_rate_prototype = zeros(2, 4))
```

In our `g!` we define the functions for computing the values of the matrix.
We can now think of the SDE that this solves as the system of equations

```math
du_1 = f_1(u,p,t)dt + g_{11}(u,p,t)dW_1 + g_{12}(u,p,t)dW_2 + g_{13}(u,p,t)dW_3 + g_{14}(u,p,t)dW_4 \\
du_2 = f_2(u,p,t)dt + g_{21}(u,p,t)dW_1 + g_{22}(u,p,t)dW_2 + g_{23}(u,p,t)dW_3 + g_{24}(u,p,t)dW_4
```

meaning that for example `du[1,1]` and `du[2,1]` correspond to stochastic changes with
the same random number in the first and second SDEs.

!!! note
    
    This problem can only be solved my SDE methods which are compatible with non-diagonal
    noise. This is discussed [in the SDE solvers page](@ref sde_solve).

The matrix itself is determined by the keyword argument `noise_rate_prototype` in the `SDEProblem`
constructor. This is a prototype for the type that `du` will be in `g!`. This can
be any `AbstractMatrix` type. Thus, we can define the problem as

```@example sde4
# Define a sparse matrix by making a dense matrix and setting some values as not zero
import SparseArrays
A = zeros(2, 4)
A[1, 1] = 1
A[1, 4] = 1
A[2, 4] = 1
A = SparseArrays.sparse(A)

# Make `g!` write the sparse matrix values
function g!(du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 4] = 0.12u[2]
    du[2, 4] = 1.8u[2]
end

# Make `g!` use the sparse matrix
prob = DE.SDEProblem(f!, g!, ones(2), (0.0, 1.0), noise_rate_prototype = A)
```

and now `g!(u,p,t)` writes into a sparse matrix, and `g!(u,p,t)*dW` is sparse matrix
multiplication.

## Example 4: Colored Noise

Colored noise can be defined [using the Noise Process interface](@ref noise_process).
In that portion of the docs, it is shown how to define your own noise process
`my_noise`, which can be passed to the SDEProblem

```julia
DE.SDEProblem(f!, g!, u0, tspan, noise = my_noise)
```

Note that general colored noise problems are only compatible with the `EM` and `EulerHeun` methods.
This is discussed [in the SDE solvers page](@ref sde_solve).

### Example: Spatially-Colored Noise in the Heston Model

Let's define the Heston equation from financial mathematics:

```math
dS = μSdt + \sqrt{v}SdW_1 \\
dv = κ(Θ-v)dt + σ\sqrt{v}dW_2 \\
dW_1 dW_2 = ρ dt
```

In this problem, we have a diagonal noise problem given by:

```@example sde4
function f!(du, u, p, t)
    du[1] = μ * u[1]
    du[2] = κ * (Θ - u[2])
end
function g!(du, u, p, t)
    du[1] = √u[2] * u[1]
    du[2] = σ * √u[2]
end
```

However, our noise has a correlation matrix for some constant `ρ`. Choosing `ρ=0.2`:

```@example sde4
ρ = 0.2
Γ = [1 ρ; ρ 1]
```

To solve this, we can define a `CorrelatedWienerProcess` which starts at zero (`W(0)=0`)
via:

```@example sde4
tspan = (0.0, 1.0)
heston_noise = DE.CorrelatedWienerProcess!(Γ, tspan[1], zeros(2), zeros(2))
```

This is then used to build the SDE:

```@example sde4
DE.SDEProblem(f!, g!, ones(2), tspan, noise = heston_noise)
```

Of course, to fully define this problem, we need to define our constants. Constructors
for making common models like this easier to define can be found in the modeling
toolkits. For example, the `HestonProblem` is pre-defined as part of the
[financial modeling tools](https://github.com/SciML/DiffEqFinancial.jl).
