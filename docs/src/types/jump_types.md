# Jump Problems

### Mathematical Specification of an problem with jumps

Jumps are defined as a Poisson process which occur according to some `rate`. When
multiple jumps are together, the process is a compound Poisson process. On their
own, a jump equation on is continuous-time Markov Chain where the time to the
next jump is exponentially distributed as calculated by the rate. This type of
process, known in biology as "Gillespie discrete stochastic simulations" and
modeled by the Chemical Master Equation (CME), is the same thing as adding jumps
to a `DiscreteProblem`. However, any differential equation can be extended by jumps
as well. For example, we have an ODE with jumps, denoted by

```math
\frac{du}{dt} = f(u,p,t) + Σ h_i(u,p,t)N_i(t)
```

where ``N_i`` is a Poisson counter of rate ``\lambda_i(u,p,t)``. Extending a stochastic
differential equation to have jumps is commonly known as a Jump Diffusion, and
is denoted by

```math
\frac{du}{dt} = f(u,p,t) + Σgᵢ(u,t)dWⁱ + Σ h_i(u,p,t)N_i(t)
```

## Variable and Constant Rate Jumps

We denote a jump as variable rate if its rate function is dependent on values which
may change between constant rate jumps. For example, if there are multiple jumps
whose rates only change when one of them occur, than that set of jumps is a constant
rate jump. If the jump's rate depends on the differential equation, time, or
by some value which changes outside of some constant rate jump, then it is denoted
as variable.

#### Defining a Constant Rate Jump

The constructor for a `ConstantRateJump` is:

```julia
ConstantRateJump(rate,affect!;save_positions=(true,true))
```

- `rate(u,p,t)` is a function which calculates the rate given the time and the state.
- `affect!(integrator)` is the effect on the equation, using the integrator interface.

#### Defining a Variable Rate Jump

The constructor for a `VariableRateJump` is:

```julia
VariableRateJump(rate,affect!;
                   idxs = nothing,
                   rootfind=true,
                   save_positions=(true,true),
                   interp_points=10,
                   abstol=1e-12,reltol=0)
```

Note that this is the same as defining a `ContinuousCallback`, except that instead
of the `condition` function, you provide a `rate(u,p,t)` function for the `rate` at
a given time and state.

## Defining a Jump Problem

To define a `JumpProblem`, you must first define the basic problem. This can be
a `DiscreteProblem` if there is no differential equation, or an ODE/SDE/DDE/DAE
if you would like to augment a differential equation with jumps. Denote this
previously defined problem as `prob`. Then the constructor for the jump problem is:

```julia
JumpProblem(prob,aggregator::Direct,jumps::JumpSet;
            save_positions = typeof(prob) <: AbstractDiscreteProblem ? (false,true) : (true,true))
```

The aggregator is the method for aggregating the constant jumps. These are defined
below. `jumps` is a `JumpSet` which is just a gathering of jumps. Instead of
passing a `JumpSet`, one may just pass a list of jumps themselves. For example:

```julia
JumpProblem(prob,aggregator,jump1,jump2)
```

and the internals will automatically build the `JumpSet`. `save_positions` is the
`save_positions` argument built by the aggregation of the constant rate jumps.

## Constant Rate Jump Aggregator

The constant rate jump aggregator is the method by which the constant rate jumps
are lumped together. This is required in all algorithms for both speed and accuracy.
The current methods are:

- `Direct`: the Gillespie SSA Direct method.

To pass the aggregator, pass the instantiation of the type. For example:

```julia
JumpProblem(prob,Direct(),jump1,jump2)
```

will build a problem where the constant rate jumps are solved using Gillespie's
Direct SSA method.
