# [Getting Started with Differential Equations in Julia](@id ode_example)

This tutorial will introduce you to the functionality for solving ODEs.
Additionally, a [video tutorial](https://youtu.be/KPEqYtEd-zY) walks through
this material.

## Example 1 : Solving Scalar Equations

In this example, we will solve the equation

```math
\frac{du}{dt} = f(u,p,t)
```

on the time interval ``t\in[0,1]`` where ``f(u,p,t)=αu``. Here, ``u`` is the
current state variable, ``p`` is our parameter variable (containing things like
a reaction rate or the constant of gravity), and ``t`` is the current time.

(In our example, we know by calculus that the solution to this equation is
``u(t)=u₀\exp(αt)``, but we will use DifferentialEquations.jl to solve this
problem *numerically*, which is essential for problems where a symbolic solution
is not known.)

The general workflow is to define a problem, solve the problem, and then analyze
the solution. The full code for solving this problem is:

```@example ODE1
import DifferentialEquations as DE
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = DE.ODEProblem(f, u0, tspan)
sol = DE.solve(prob, DE.Tsit5(), reltol = 1e-8, abstol = 1e-8)

import Plots
Plots.plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
Plots.plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")
```

where the pieces are described below.

### Step 1: Defining a Problem

To solve this numerically, we define a problem type by giving it the equation,
the initial condition, and the timespan to solve over:

```@example ODE2
import DifferentialEquations as DE
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = DE.ODEProblem(f, u0, tspan)
```

Note that DifferentialEquations.jl will choose the types for the problem based on
the types used to define the problem type. For our example, notice that `u0` is
a Float64, and therefore this will solve with the dependent variables being
Float64. Since `tspan = (0.0,1.0)` is a tuple of Float64's, the independent variables
will be solved using Float64's (note that the start time and end time must match
types). You can use this to choose to solve with arbitrary precision numbers,
unitful numbers, etc. Please see the
[notebook tutorials](https://github.com/SciML/SciMLTutorials.jl)
for more examples.

The problem types include many other features, including the ability to define
mass matrices and hold callbacks for events. Each problem type has a page which
details its constructor and the available fields. For
[ODEs, the appropriate page is here](@ref ode_prob).
In addition, a user can specify additional functions to be associated with the
function in order to speed up the solvers. These are detailed
[at the performance overloads page](@ref performance_overloads).

### Step 2: Solving a Problem

#### Controlling the Solvers

After defining a problem, you solve it using `solve`.

```@example ODE2
sol = DE.solve(prob)
```

This gives us an object `sol` which contains the solution. Looking at the solution object:

```@example ODE2
typeof(sol)
```

The solution object contains the time points and corresponding solution values:

```@example ODE2
@info "Solution contains $(length(sol.t)) time points from t=$(sol.t[1]) to t=$(sol.t[end])"
```

The solvers can be controlled using the available options are described on the
[Common Solver Options manual page](@ref solver_options). For example,
we can lower the relative tolerance (in order to get a more correct result, at
the cost of more timesteps) by using the command `reltol`:

```@example ODE2
sol = DE.solve(prob, reltol = 1e-6);
nothing # hide
```

There are many controls for handling outputs. For example, we can choose to have
the solver save every `0.1` time points by setting `saveat=0.1`. Chaining this
with the tolerance choice looks like:

```@example ODE2
sol = DE.solve(prob, reltol = 1e-6, saveat = 0.1);
nothing # hide
```

More generally, `saveat` can be any collection of time points to save at.
Note that this uses interpolations to keep the timestep unconstrained to speed
up the solution. In addition, if we only care about the endpoint, we can turn
off intermediate saving in general:

```@example ODE2
sol = DE.solve(prob, reltol = 1e-6, save_everystep = false);
nothing # hide
```

which will only save the final time point.

#### Choosing a Solver Algorithm

DifferentialEquations.jl has a method for choosing the default solver algorithm,
which will find an efficient method to solve your problem. To help users receive
the right algorithm, DifferentialEquations.jl offers a method for choosing
algorithms through hints. This default chooser utilizes the precision of the
number types and the keyword arguments (such as the tolerances) to select an
algorithm. Additionally, one can provide `alg_hints` to help choose good defaults
using properties of the problem and necessary features for the solution.
For example, if we have a stiff problem where we need high accuracy,
but don't know the best stiff algorithm for this problem, we can use:

```@example ODE2
sol = DE.solve(prob, alg_hints = [:stiff], reltol = 1e-8, abstol = 1e-8);
nothing # hide
```

You can also explicitly choose the algorithm to use. DifferentialEquations.jl
offers a much wider variety of solver algorithms than traditional differential
equations libraries. Many of these algorithms are from recent research and have
been shown to be more efficient than the “standard” algorithms.
For example, we can choose a 5th order Tsitouras method:

```@example ODE2
sol = DE.solve(prob, DE.Tsit5());
nothing # hide
```

Note that the solver controls can be combined with the algorithm choice. Thus
we can for example solve the problem using `DE.Tsit5()` with a lower tolerance
via:

```@example ODE2
sol = DE.solve(prob, DE.Tsit5(), reltol = 1e-8, abstol = 1e-8);
nothing # hide
```

In DifferentialEquations.jl, some good “go-to” choices for ODEs are:

  - `AutoTsit5(Rosenbrock23())` handles both stiff and non-stiff equations. This
    is a good algorithm to use if you know nothing about the equation.
  - `AutoVern7(Rodas5())` handles both stiff and non-stiff equations in a way that's
    efficient for high accuracy.
  - `DE.Tsit5()` for standard non-stiff. This is the first algorithm to try in
    most cases.
  - `BS3()` for fast low accuracy non-stiff.
  - `Vern7()` for high accuracy non-stiff.
  - `Rodas4()` or `Rodas5()` for small stiff equations with Julia-defined types, events, etc.
  - `KenCarp4()` or `TRBDF2()` for medium-sized (100-2000 ODEs) stiff equations
  - `RadauIIA5()` for really high accuracy stiff equations
  - `QNDF()` for large stiff equations

For a comprehensive list of the available algorithms and detailed recommendations,
[please see the solver documentation](@ref ode_solve). Every problem
type has an associated page detailing all the solvers associated with the problem.

### Step 3: Analyzing the Solution

#### Handling the Solution Type

The result of `solve` is a solution object. We can access the 5th value of the
solution with:

```@example ODE2
sol[5]
```

or get the time of the 8th timestep by:

```@example ODE2
sol.t[8]
```

Convenience features are also included. We can build an array using a
comprehension over the solution tuples via:

```@example ODE2
[t + u for (u, t) in DE.tuples(sol)]
```

or more generally

```@example ODE2
[t + 2u for (u, t) in zip(sol.u, sol.t)]
```

allows one to use more parts of the solution type. The object that is returned by
default acts as a continuous solution via an interpolation. We can access the
interpolated values by treating `sol` as a function, for example:

```@example ODE2
sol(0.45) # The value of the solution at t=0.45
```

Note the difference between these: indexing with `[i]` is the value at the `i`th
step, while `(t)` is an interpolation at time `t`!

If in the solver `dense=true` (this is the default unless `saveat` is used), then
this interpolation is a high order interpolation and thus usually matches the
error of the solution time points. The interpolations associated with each solver
is [detailed at the solver algorithm page](@ref ode_solve). If `dense=false`
(unless specifically set, this only occurs when `save_everystep=false` or `saveat`
is used) then this defaults to giving a linear interpolation.

For more details on handling the output, see [the solution handling page](@ref solution).

#### Plotting Solutions

While one can directly plot solution time points using the tools given above,
convenience commands are defined by recipes for Plots.jl. To plot the solution
object, simply call plot:

```@example ODE2
#]add Plots # You need to install Plots.jl before your first time using it!
import Plots
#plotly() # You can optionally choose a plotting backend
Plots.plot(sol)
```

The plot function can be formatted using [the attributes available in Plots.jl](https://juliaplots.org/).
Additional DiffEq-specific controls are documented [at the plotting page](@ref plot).

For example, from the Plots.jl attribute page, we see that the line width can be
set via the argument `linewidth`. Additionally, a title can be set with `title`.
Thus we add these to our plot command to get the correct output, fix up some
axis labels, and change the legend (note we can disable the legend with
`legend=false`) to get a nice-looking plot:

```@example ODE2
Plots.plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
```

We can then add to the plot using the `plot!` command:

```@example ODE2
Plots.plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")
```

## Example 2: Solving Systems of Equations

In this example, we will solve the Lorenz equations:

```math
\begin{aligned}
\frac{dx}{dt} &= σ(y-x) \\
\frac{dy}{dt} &= x(ρ-z) - y \\
\frac{dz}{dt} &= xy - βz \\
\end{aligned}
```

Defining your ODE function to be in-place updating can have performance benefits.
What this means is that, instead of writing a function which outputs its solution,
you write a function which updates a vector that is designated to hold the solution.
By doing this, DifferentialEquations.jl's solver packages are able to reduce the
amount of array allocations and achieve better performance.

The way we do this is we simply write the output to the 1st input of the function.
For example, our Lorenz equation problem would be defined by the function:

```@example ODE3
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
```

and then we can use this function in a problem:

```@example ODE3
import DifferentialEquations as DE
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = DE.ODEProblem(lorenz!, u0, tspan)
sol = DE.solve(prob)
@info "Solution has $(length(sol.t)) timesteps" # hide
nothing # hide
```

Using the plot recipe tools
[defined on the plotting page](@ref plot_vars),
we can choose to do a 3D phase space plot between the different variables:

```@example ODE3
import Plots
Plots.plot(sol, idxs = (1, 2, 3))
```

Note that the default plot for multidimensional systems is an overlay of
each timeseries. We can plot the timeseries of just the second component using
the variable choices interface once more:

```@example ODE3
Plots.plot(sol, idxs = (0, 2))
```

Note that here “variable 0” corresponds to the independent variable (“time”).

## Defining Parameterized Functions

Often, you want to explicitly have parameters associated with your
differential equations. This can be used by things like
[parameter estimation routines](https://docs.sciml.ai/Overview/stable/highlevels/inverse_problems/).
In this case, you use the `p` values via the syntax:

```@example ODE3
function parameterized_lorenz!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
end
```

and then we add the parameters to the `ODEProblem`:

```@example ODE3
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 1.0)
p = [10.0, 28.0, 8 / 3]
prob = DE.ODEProblem(parameterized_lorenz!, u0, tspan, p)
```

We can make our functions look nicer by doing a few tricks. For example:

```@example ODE3
function parameterized_lorenz!(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = dx = σ * (y - x)
    du[2] = dy = x * (ρ - z) - y
    du[3] = dz = x * y - β * z
end
```

Note that the type for the parameters `p` can be anything: you can use arrays,
static arrays, named tuples, etc. to enclose your parameters in a way that is
sensible for your problem.

Since the parameters exist within the function, functions defined in this manner
can also be used for sensitivity analysis, parameter estimation routines,
and bifurcation plotting. This makes DifferentialEquations.jl a full-stop solution
for differential equation analysis which also achieves high performance.

## Example 3: Solving Nonhomogeneous Equations using Parameterized Functions

Parameterized functions can also be used for building **nonhomogeneous ordinary differential equations** (these are also referred to as ODEs with **nonzero right-hand sides**). They are frequently used as models for dynamical systems with external (in general time-varying) **inputs**. As an example, consider a [model of a pendulum](https://en.wikipedia.org/wiki/Pendulum_(mathematics)) consisting of a slender rod of length `l` and mass `m`:

```math
\begin{aligned}
\frac{\mathrm{d}\theta(t)}{\mathrm{d}t} &= \omega(t)\\
\frac{\mathrm{d}\omega(t)}{\mathrm{d}t} &= - \frac{3}{2}\frac{g}{l}\sin\theta(t) + \frac{3}{ml^2}M(t)
\end{aligned},
```

where `θ` and `ω` are the angular deviation of the pendulum from the vertical (hanging) orientation and the angular rate, respectively, `M` is an external torque (developed, say, by a wind or a motor), and finally, `g` stands for gravitational acceleration.

```@example ODE4
import DifferentialEquations as DE
import Plots

l = 1.0                             # length [m]
m = 1.0                             # mass [kg]
g = 9.81                            # gravitational acceleration [m/s²]

function pendulum!(du, u, p, t)
    du[1] = u[2]                    # θ'(t) = ω(t)
    du[2] = -3g / (2l) * sin(u[1]) + 3 / (m * l^2) * p(t) # ω'(t) = -3g/(2l) sin θ(t) + 3/(ml^2)M(t)
end

θ₀ = 0.01                           # initial angular deflection [rad]
ω₀ = 0.0                            # initial angular velocity [rad/s]
u₀ = [θ₀, ω₀]                       # initial state vector
tspan = (0.0, 10.0)                  # time interval

M = t -> 0.1sin(t)                    # external torque [Nm]

prob = DE.ODEProblem(pendulum!, u₀, tspan, M)
sol = DE.solve(prob)

Plots.plot(sol, linewidth = 2, xaxis = "t", label = ["θ [rad]" "ω [rad/s]"], layout = (2, 1))
```

Note how the external **time-varying** torque `M` is introduced as a **parameter** in the `pendulum!` function. Indeed, as a general principle the parameters can be any type; here we specify `M` as time-varying by representing it by a function, which is expressed by appending the dependence on time `(t)` to the name of the parameter.

Note also that, in contrast with the time-varying parameter, the (vector of) state variables `u`, which is generally also time-varying, is always used without the explicit dependence on time `(t)`.

## [Example 4: Using Other Types for Systems of Equations](@id ode_other_types)

DifferentialEquations.jl can handle many different dependent variable types
(generally, anything with a linear index should work!). So instead of solving a
vector equation, let's let `u` be a matrix! To do this, we simply need to have `u0`
be a matrix, and define `f` such that it takes in a matrix and outputs a matrix.
We can define a matrix of linear ODEs as follows:

```@example ODE4
import DifferentialEquations as DE
import Plots
A = [1.0 0 0 -5
     4 -2 4 -3
     -4 0 0 1
     5 -2 2 3]
u0 = rand(4, 2)
tspan = (0.0, 1.0)
f(u, p, t) = A * u
prob = DE.ODEProblem(f, u0, tspan)
```

Here our ODE is on a 4x2 matrix, and the ODE is the linear system defined by
multiplication by `A`. To solve the ODE, we do the same steps
as before.

```@example ODE4
sol = DE.solve(prob)
Plots.plot(sol)
```

We can instead use the in-place form by using Julia's in-place matrix multiplication
function `mul!`:

```@example ODE4
import LinearAlgebra
f(du, u, p, t) = LinearAlgebra.mul!(du, A, u)
```

Additionally, we can use non-traditional array types as well. For example,
`StaticArrays.jl` offers immutable arrays which are stack-allocated, meaning
that their usage does not require any (slow) heap-allocations that arrays
normally have. This means that they can be used to solve the same problem as
above, with the only change being the type for the initial condition and constants:

```@example ODE4
import StaticArrays
A = StaticArrays.@SMatrix [1.0 0.0 0.0 -5.0
              4.0 -2.0 4.0 -3.0
              -4.0 0.0 0.0 1.0
              5.0 -2.0 2.0 3.0]
u0 = StaticArrays.@SMatrix rand(4, 2)
tspan = (0.0, 1.0)
f2(u, p, t) = A * u
prob = DE.ODEProblem(f2, u0, tspan)
sol = DE.solve(prob)
Plots.plot(sol)
```

Note that the analysis tools generalize over to systems of equations as well.

```@example ODE4
sol[4]
```

still returns the solution at the fourth timestep. It also indexes into the array
as well. The last value is the timestep, and the beginning values are for the component.
This means

```@example ODE4
sol[5, 3]
```

is the value of the 5th component (by linear indexing) at the 3rd timepoint, or

```@example ODE4
sol[2, 1, :]
```

is the timeseries for the component, which is the 2nd row and 1 column.

## Going Beyond ODEs: How to Use the Documentation

Not everything can be covered in the tutorials. Instead, this tutorial will end
by pointing you in the directions for the next steps.

#### Common API for Defining, Solving, and Plotting

One feature of DifferentialEquations.jl is that this pattern for solving equations
is conserved across the different types of differential equations. Every equation
has a problem type, a solution type, and the same solution handling (+ plotting)
setup. Thus the solver and plotting commands in the **Basics** section applies to
all sorts of equations, like stochastic differential equations and delay differential
equations. Each of these different problem types are defined in the **Problem Types**
section of the docs. Every associated solver algorithm is detailed in the
**Solver Algorithms** section, sorted by problem type. The same steps for ODEs
can then be used for the analysis of the solution.

#### Additional Features and Analysis Tools

In many cases, the common workflow only starts with solving the differential equation.
Many common setups have built-in solutions in DifferentialEquations.jl. For example,
check out the features for:

  - [Handling, parallelizing, and analyzing large Ensemble experiments](@ref ensemble)
  - [Saving the output to tabular formats like DataFrames and CSVs](@ref io)
  - [Event handling](@ref callbacks)
  - [Parameter estimation (inverse problems)](https://docs.sciml.ai/Overview/stable/highlevels/inverse_problems/)
  - [Quantification of numerical uncertainty and error](https://docs.sciml.ai/Overview/stable/highlevels/uncertainty_quantification/)

Many more are defined in the relevant sections of the docs. Please explore the rest
of the documentation, including tutorials for getting started with other types
of equations. In addition, to get help, please either
[file an issue at the main repository](https://github.com/SciML/DifferentialEquations.jl)
or [come have an informal discussion at the Julia Zulip chatroom](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged).
