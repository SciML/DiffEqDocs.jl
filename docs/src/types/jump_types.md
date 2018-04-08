# Jump Problems

### Mathematical Specification of an problem with jumps

Jumps are defined as a Poisson process which changes states at some `rate`. When
there are multiple possible jumps, the process is a compound Poisson process. On their
own, a jump equation is continuous-time Markov Chain where the time to the
next jump is exponentially distributed as calculated by the rate. This type of
process, known in biology as "Gillespie discrete stochastic simulations" and
modeled by the Chemical Master Equation (CME), is the same thing as adding jumps
to a `DiscreteProblem`. However, any differential equation can be extended by jumps
as well. For example, we have an ODE with jumps, denoted by

```math
\frac{du}{dt} = f(u,p,t) + Σ c_i(u,p,t)dp_i
```

where ``dp_i`` is a Poisson counter of rate ``\lambda_i(u,p,t)``. Extending a
stochastic differential equation to have jumps is commonly known as a Jump
Diffusion, and is denoted by

```math
\frac{du}{dt} = f(u,p,t) + Σgᵢ(u,t)dWⁱ + Σ c_i(u,p,t)dp_i
```

## Types of Jumps: Regular, Variable, Constant Rate and Mass Action

A `RegularJump` is a set of jumps that do not make structural changes to the
underlying equation. These kinds of jumps only change values of the dependent
variable (`u`) and thus can be treated in an inexact manner. Other jumps, such
as those which change the size of `u`, require exact handling which is also
known as time-adaptive jumping. These can only be specified as a
`ConstantRateJump`, `MassActionJump`, or a `VariableRateJump`.

We denote a jump as variable rate if its rate function is dependent on values
which may change between constant rate jumps. For example, if there are multiple
jumps whose rates only change when one of them occur, than that set of jumps is
a constant rate jump. If a jump's rate depends on the differential equation,
time, or by some value which changes outside of any constant rate jump, then it
is denoted as variable.

A `MassActionJump` is a specialized representation for a collection of constant
rate jumps that can each be interpreted as a standard mass action reaction. For
systems comprised of many mass action reactions, using the `MassActionJump` type
will offer improved performance. Note, only one `MassActionJump` should be
defined per `JumpProblem`; it is then responsible for handling all mass action
reaction type jumps. For systems with both mass action jumps and non-mass action
jumps, one can create one `MassActionJump` to handle the mass action jumps, and
create a number of `ConstantRateJumps` to handle the non-mass action jumps.

`RegularJump`s are optimized for regular jumping algorithms like tau-leaping and
hybrid algorithms. `ConstantRateJump`s and `MassActionJump`s are optimized for
SSA algorithms. `ConstantRateJump`s, `MassActionJump`s and `VariableRateJump`s
can be added to standard DiffEq algorithms since they are simply callbacks,
while `RegularJump`s require special algorithms. 

#### Defining a Regular Jump

The constructor for a `RegularJump` is:

```julia
RegularJump(rate,c,c_prototype;mark_dist = nothing,constant_c = false)
```

- `rate(out,u,p,t)` is the function which computes the rate for every regular
  jump process
- `c(dc,u,p,t,mark)` is the current Stoichiometry matrix for each jump process
- `dc` is the cache array to be used for `dc`
- `mark_dist` is the distribution for the mark
- `constant_c` denotes whether the Stoichiometry matrix `c` is constant

`dc` is an `n x m` matrix, where `n` is the number of Poisson processes and `m`
is the number of dependent variables (should match `length(u)`). `rate` is a
vector equation which should compute the rates in to `out` which is a length
`n` vector.

#### Defining a Constant Rate Jump

The constructor for a `ConstantRateJump` is:

```julia
ConstantRateJump(rate,affect!)
```

- `rate(u,p,t)` is a function which calculates the rate given the time and the state.
- `affect!(integrator)` is the effect on the equation, using the integrator interface.


#### Defining a Mass Action Jump

The constructor for a `MassActionJump` is:
```julia
MassActionJump(rate_consts, reactant_stoich, net_stoich; scale_rates = true)
```
- `rate_consts` is a vector of the rate constants for each reaction.
- `reactant_stoich` is a vector whose `k`th entry is the reactant stoichiometry
  of the `k`th reaction. The reactant stoichiometry for an individual reaction
  is assumed to be represented as a vector of `Pair`s, mapping species id to
  stoichiometric coefficient.
- `net_stoich` is assumed to have the same type as `reactant_stoich`; a
  vector whose `k`th entry is the net stoichiometry of the `k`th reaction. The
  net stoichiometry for an individual reaction is again represented as a vector
  of `Pair`s, mapping species id to the net change in the species when the
  reaction occurs.
- `scale_rates` is an optional parameter that specifies whether the rate
  constants correspond to stochastic rate constants in the sense used by
  Gillespie, and hence need to be rescaled. *The default, `scale_rates=true`,
  corresponds to rescaling the passed in rate constants.* See below.

When using `MassActionJump` the default behavior is to assume rate constants
correspond to stochastic rate constants in the sense used by Gillespie (J. Comp.
Phys., 1976, 22 (4)). This means that for a reaction such as
``2A \overset{k}{\rightarrow} B``, the jump rate function constructed by
`MassActionJump` would be `k*A*(A-1)/2!`. For a trimolecular reaction like ``3A
\overset{k}{\rightarrow} B`` the rate function would be `k*A*(A-1)*(A-2)/3!`. To *avoid*
having the reaction rates rescaled (by `1/2` and `1/6` for these two examples),
one can pass the `MassActionJump` constructor the optional named parameter
`scale_rates=false`, i.e. use
```julia
MassActionJump(rates, reactant_stoich, net_stoich; scale_rates = false)
```


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

Note that a `JumpProblem`/`JumpSet` can only have 1 `RegularJump` (since a
`RegularJump` itself describes multiple processes together). Similarly, it can
only have one `MassActionJump` (since it also describes multiple processes
together).

## Constant Rate Jump Aggregators

Constant rate jump aggregators are the methods by which constant rate
jumps, including `MassActionJump`s, are lumped together. This is required in all
algorithms for both speed and accuracy. The current methods are:

- `Direct`: the Gillespie Direct method SSA.
- `DirectFW`: the Gillespie Direct method SSA with `FunctionWrappers`. This
  aggregator uses a different internal storage format for collections of
  `ConstantRateJumps`. For systems with ~10 or more `ConstantRateJumps` it
  should offer better performance than `Direct`.
- `FRM`: the Gillespie first reaction method SSA. `Direct` should generally
  offer better performance and be preferred to `FRM`.
- `FRMFW`: the Gillespie first reaction method SSA with `FunctionWrappers`.
  `DirectFW` should generally offer better performance and be preferred to `FRMFW`.


To pass the aggregator, pass the instantiation of the type. For example:

```julia
JumpProblem(prob,Direct(),jump1,jump2)
```

will build a problem where the constant rate jumps are solved using Gillespie's
Direct SSA method.


## Recommendations for Constant Rate Jumps
For representing and aggregating constant rate jumps 
- Use a `MassActionJump` to handle all jumps that can be represented as mass
  action reactions. This will generally offer the fastest performance. 
- Use `ConstantRateJump`s for any remaining jumps.
  - If there are *less* than ~10 `ConstantRateJumps`, the `Direct` aggregator
    will generally offer the best performance.
  - If there are *more* than ~10 `ConstantRateJumps`, the `DirectFW` aggregator
    will generally offer the best performance.
