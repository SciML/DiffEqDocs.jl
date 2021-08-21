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
RegularJump(rate,c,numjumps;mark_dist = nothing)
```

- `rate(out,u,p,t)` is the function which computes the rate for every regular
  jump process
- `c(du,u,p,t,counts,mark)` is calculates the update given `counts` number of
  jumps for each jump process in the interval.
- `numjumps` is the number of jump processes, i.e. the number of `rate` equations
  and the number of `counts`
- `mark_dist` is the distribution for the mark.

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
MassActionJump(reactant_stoich, net_stoich; scale_rates = true, param_idxs=nothing)
```
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
- `param_idxs` is a vector of the indices within the parameter vector, `p`, that
  correspond to the rate constant for each jump.

**Notes for Mass Action Jumps**
- When using `MassActionJump` the default behavior is to assume rate constants
  correspond to stochastic rate constants in the sense used by Gillespie (J.
  Comp. Phys., 1976, 22 (4)). This means that for a reaction such as ``2A
  \overset{k}{\rightarrow} B``, the jump rate function constructed by
  `MassActionJump` would be `k*A*(A-1)/2!`. For a trimolecular reaction like
  ``3A \overset{k}{\rightarrow} B`` the rate function would be
  `k*A*(A-1)*(A-2)/3!`. To *avoid* having the reaction rates rescaled (by `1/2`
  and `1/6` for these two examples), one can pass the `MassActionJump`
  constructor the optional named parameter `scale_rates=false`, i.e. use
  ```julia
  MassActionJump(reactant_stoich, net_stoich; scale_rates = false, param_idxs)
  ```
- Zero order reactions can be passed as `reactant_stoich`s in one of two ways.
  Consider the ``\varnothing \overset{k}{\rightarrow} A`` reaction with rate `k=1`:
  ```julia
  p = [1.]
  reactant_stoich = [[0 => 1]]
  net_stoich = [[1 => 1]]
  jump = MassActionJump(reactant_stoich, net_stoich; param_idxs=[1])
  ```
  Alternatively one can create an empty vector of pairs to represent the reaction:
  ```julia
  p = [1.]
  reactant_stoich = [Vector{Pair{Int,Int}}()]
  net_stoich = [[1 => 1]]
  jump = MassActionJump(reactant_stoich, net_stoich; param_idxs=[1])
  ```
- For performance reasons, it is recommended to order species indices in
  stoichiometry vectors from smallest to largest. That is 
  ```julia
  reactant_stoich = [[1 => 2, 3 => 1, 4 => 2], [2 => 2, 3 => 2]]
  ```
  is preferred over
  ```julia
  reactant_stoich = [[3 => 1, 1 => 2, 4 = > 2], [3 => 2, 2 => 2]]
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
- `RDirect`: A variant of Gillespie's Direct method that uses rejection to
  sample the next reaction.
- *`DirectCR`*: The Composition-Rejection Direct method of Slepoy et al. For
  large networks and linear chain-type networks it will often give better
  performance than `Direct`. (Requires dependency graph, see below.)
- `DirectFW`: the Gillespie Direct method SSA with `FunctionWrappers`. This
  aggregator uses a different internal storage format for collections of
  `ConstantRateJumps`. 
- `FRM`: the Gillespie first reaction method SSA. `Direct` should generally
  offer better performance and be preferred to `FRM`.
- `FRMFW`: the Gillespie first reaction method SSA with `FunctionWrappers`.
- *`NRM`*: The Gibson-Bruck Next Reaction Method. For some reaction network
   structures this may offer better performance than `Direct` (for example,
   large, linear chains of reactions). (Requires dependency graph, see below.) 
- *`RSSA`*: The Rejection SSA (RSSA) method of Thanh et al. With `RSSACR`, for
  very large reaction networks it often offers the best performance of all
  methods. (Requires dependency graph, see below.)
- *`RSSACR`*: The Rejection SSA (RSSA) with Composition-Rejection method of
  Thanh et al. With `RSSA`, for very large reaction networks it often offers the
  best performance of all methods. (Requires dependency graph, see below.)
- *`SortingDirect`*: The Sorting Direct Method of McCollum et al. It will
  usually offer performance as good as `Direct`, and for some systems can offer
  substantially better performance. (Requires dependency graph, see below.)

To pass the aggregator, pass the instantiation of the type. For example:

```julia
JumpProblem(prob,Direct(),jump1,jump2)
```

will build a problem where the constant rate jumps are solved using Gillespie's
Direct SSA method.

## Constant Rate Jump Aggregators Requiring Dependency Graphs
Italicized constant rate jump aggregators require the user to pass a dependency
graph to `JumpProblem`. `DirectCR`, `NRM` and `SortingDirect` require a
jump-jump dependency graph, passed through the named parameter `dep_graph`. i.e.
```julia
JumpProblem(prob,DirectCR(),jump1,jump2; dep_graph=your_dependency_graph)
```
For systems with only `MassActionJump`s, or those generated from a
[Catalyst](https://github.com/SciML/Catalyst.jl) `reaction_network`, this graph
will be auto-generated. Otherwise you must construct the dependency graph
manually. Dependency graphs are represented as a `Vector{Vector{Int}}`, with the
`i`th vector containing the indices of the jumps for which rates must be
recalculated when the `i`th jump occurs. Internally, all `MassActionJump`s are
ordered before `ConstantRateJump`s (with the latter internally ordered in the
same order they were passed in).

`RSSA` and `RSSACR` require two different types of dependency graphs, passed
through the following `JumpProblem` kwargs:
- `vartojumps_map` - A `Vector{Vector{Int}}` mapping each variable index, `i`,
  to a set of jump indices. The jump indices correspond to jumps with rate
  functions that depend on the value of `u[i]`.
-  `jumptovars_map` - A `Vector{Vector{Int}}`  mapping each jump index to a set
   of variable indices. The corresponding variables are those that have their
   value, `u[i]`, altered when the jump occurs.

For systems generated from a [Catalyst](https://github.com/SciML/Catalyst.jl)
`reaction_network` these will be auto-generated. Otherwise you must explicitly
construct and pass in these mappings.

## Recommendations for Constant Rate Jumps
For representing and aggregating constant rate jumps 
- Use a `MassActionJump` to handle all jumps that can be represented as mass
  action reactions. This will generally offer the fastest performance. 
- Use `ConstantRateJump`s for any remaining jumps.
- For a small number of jumps, < ~10, `Direct` will often perform as well as the
  other aggregators.
- For > ~10 jumps `SortingDirect` will often offer better performance than `Direct`.
- For large numbers of jumps with sparse chain like structures and similar jump
  rates, for example continuous time random walks, `RSSACR`, `DirectCR` and then
  `NRM` often have the best performance.
- For very large networks, with many updates per jump, `RSSA` and `RSSACR` will
  often substantially outperform the other methods. 

In general, for systems with sparse dependency graphs if `Direct` is slow, one
of `SortingDirect`, `RSSA` or `RSSACR` will usually offer substantially better
performance. See
[DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl) for
benchmarks on several example networks.

## Remaking `JumpProblem`s
When running many simulations, it can often be convenient to update the initial
condition or simulation parameters without having to create and initialize a new
`JumpProblem`. In such situations `remake` can be used to change the initial
condition, time span, and the parameter vector. **Note,** the new `JumpProblem`
will alias internal data structures from the old problem, including core
components of the SSA aggregators. As such, only the new problem generated by
`remake` should be used for subsequent simulations.

As an example, consider the following SIR model:
```julia
rate1(u,p,t) = (0.1/1000.0)*u[1]*u[2]
function affect1!(integrator)
  integrator.u[1] -= 1
  integrator.u[2] += 1
end
jump = ConstantRateJump(rate1,affect1!)

rate2(u,p,t) = 0.01u[2]
function affect2!(integrator)
  integrator.u[2] -= 1
  integrator.u[3] += 1
end
jump2 = ConstantRateJump(rate2,affect2!)
u0    = [999,1,0]
p     = (0.1/1000,0.01)
tspan = (0.0,250.0)
dprob = DiscreteProblem(u0, tspan, p)
jprob = JumpProblem(dprob, Direct(), jump, jump2)
sol   = solve(jprob, SSAStepper())
```
We can change any of `u0`, `p` and `tspan` by either making a new
`DiscreteProblem`
```julia
u02    = [10,1,0]
p2     = (.1/1000, 0.0)
tspan2 = (0.0,2500.0)
dprob2 = DiscreteProblem(u02, tspan2, p2)
jprob2 = remake(jprob, prob=dprob2)
sol2   = solve(jprob2, SSAStepper())
```
or by directly remaking with the new parameters
```julia
jprob2 = remake(jprob, u0=u02, p=p2, tspan=tspan2)
sol2   = solve(jprob2, SSAStepper())
```
To avoid ambiguities, the following will give an error
```julia
jprob2 = remake(jprob, prob=dprob2, u0=u02)
```
as will trying to update either `p` or `tspan` while passing a new
`DiscreteProblem` using the `prob` kwarg.

