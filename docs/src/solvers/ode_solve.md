# ODE Solvers

`solve(prob::ODEProblem,alg;kwargs)`

Solves the ODE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

It is suggested that you try choosing an algorithm using the `alg_hints`
keyword argument. However, in some cases you may want something specific,
or you may just be curious. This guide is to help you choose the right algorithm.

### Unknown Stiffness Problems

When the stiffness of the problem is unknown, it is recommended you use a
stiffness detection and auto-switching algorithm. These methods are multi-paradigm
and allow for efficient solution of both stiff and non-stiff problems. The cost
for auto-switching is very minimal but the choices are restrained and so they
are a good go-to method when applicable.

For default tolerances, `AutoTsit5(Rosenbrock23())` is a good choice. For lower
tolerances, using `AutoVern7` or `AutoVern9` with `Rodas4`, `KenCarp4`, or
`Rodas5` can all be good choices depending on the problem. For very large
systems (`>1000 ODEs?`), consider using `lsoda`.

### Non-Stiff Problems

For non-stiff problems, the native OrdinaryDiffEq.jl algorithms are vastly
more efficient than the other choices. For most non-stiff
problems, we recommend `Tsit5`. When more robust error control is required,
`BS5` is a good choice. If at moderate tolerances an the interpolation error
is very important, consider the `OwrenZen5` method. For fast solving at higher
tolerances, we recommend `BS3` (or `OwrenZen3` if the interpolation error is
important). For high accuracy but with the range of `Float64` (`~1e-8-1e-12`),
we recommend `Vern6`, `Vern7`, or `Vern8` as efficient choices.

For high accuracy non-stiff solving (`BigFloat` and tolerances like `<1e-12`),
we recommend the `Vern9` method. If a high-order method is needed with a high
order interpolant, then you should choose `Vern9` which is Order 9 with an
Order 9 interpolant. If you need extremely high accuracy (`<1e-30`?) and do
not need an interpolant, try the `Feagin12` or `Feagin14` methods. Note that the
Feagin methods are the only high-order optimized methods which do not include a
high-order interpolant (they do include a 3rd order Hermite interpolation if
needed). Note that these high order RK methods are more robust than the high order
Adams-Bashforth methods to discontinuities and achieve very high precision, and
are much more efficient than the extrapolation methods. However, the `CVODE_Adams`
method can be a good choice for high accuracy when the system of equations is
very large (`>1,000` ODEs?), the function calculation is very expensive,
or the solution is very smooth.

If strict error bounds are needed, then adaptive methods with defect controls
are required. Defect controls use an error measurement on the interpolating
polynomial to make the error estimate better capture the error over the full
interval. For medium accuracy calculations, `RK4` is a good choice.

### Stiff Problems

For stiff problems at high tolerances (`>1e-2`?) it is recommended that you use
`Rosenbrock23` or `TRBDF2`. These are robust to oscillations and massive
stiffness is needed, though are only efficient when low accuracy is needed.
`Rosenbrock23` is more efficient for small systems where re-evaluating and
re-factorizing the Jacobian is not too costly, and for sufficiently large
systems `TRBDF2` will be more efficient. `ABDF2` can be the most efficient
the largest systems or most expensive `f`.

At medium tolerances (`>1e-8`?) it is recommended you use `Rodas5`,
`Rodas4P` (the former is more efficient but the later is more reliable),
`Kvaerno5`, or `KenCarp4`. As native DifferentialEquations.jl solvers,
many Julia numeric types (such as BigFloats,
[ArbFloats](https://github.com/JuliaArbTypes/ArbFloats.jl), or
[DecFP](https://github.com/stevengj/DecFP.jl)) will work. When the equation is
defined via the `@ode_def` macro, these will be the most efficient.

For faster solving at low tolerances (`<1e-9`) but when `Vector{Float64}` is used,
use `radau`.

For asymptotically large systems of ODEs (`N>1000`?)
where `f` is very costly and the complex eigenvalues are minimal (low oscillations),
in that case `CVODE_BDF` will be the most efficient but requires `Vector{Float64}`.
`CVODE_BDF` will also do surprisingly well if the solution is smooth. However,
this method can be less stiff than other methods and stuff may fail at low
accuracy situations. Another good choice for this regime is `lsoda`.

#### Special Properties of Stiff Integrators

`ImplicitMidpoint` is a symmetric and symplectic integrator.
`Trapezoid` is a symmetric (almost symplectic) integrator with adaptive
timestepping. `ImplicitEuler` is an extension to the common algorithm with
adaptive timestepping and efficient quasi-Newton Jacobian re-usage which is fully
strong-stability presurving (SSP) for hyperbolic PDEs.

Notice that `Rodas4` loses accuracy on discretizations of nonlinear
parabolic PDEs, and thus it's suggested you replace it with `Rodas4P` in those
situations which is 3rd order. `ROS3P` is only third order and achieves 3rd order
on such problems and can thus be more efficient in this case.

## Translations from MATLAB/Python/R

For users familiar with MATLAB/Python/R, good translations of the standard
library methods are as follows:

- `ode23` --> `BS3()`
- `ode45`/`dopri5` --> `DP5()`, though in most cases `Tsit5()` is more efficient
- `ode23s` --> `Rosenbrock23()`, though in most cases `Rodas4()` is more efficient
- `ode113` --> `CVODE_Adams()`, though in many cases `Vern7()` is more efficient
- `dop853` --> `DP8()`, though in most cases `Vern7()` is more efficient
- `ode15s`/`vode` --> `CVODE_BDF()`, though in many cases `Rodas4()` or `radau()`
  are more efficient
- `ode23t` --> `Trapezoid()` for efficiency and `GenericTrapezoid()` for robustness
- `ode23tb` --> `TRBDF2`
- `lsoda` --> `lsoda()` (requires `Pkg.add("LSODA"); using LSODA`)
- `ode15i` --> `IDA()`, though in many cases `Rodas4()` can handle the DAE and is
  significantly more efficient

# Full List of Methods

Choose one of these methods with the `alg` keyword in `solve`.

## OrdinaryDiffEq.jl

Unless otherwise specified, the OrdinaryDiffEq algorithms all come with a
3rd order Hermite polynomial interpolation. The algorithms denoted as having a "free"
interpolation means that no extra steps are required for the interpolation. For
the non-free higher order interpolating functions, the extra steps are computed
lazily (i.e. not during the solve).

The OrdinaryDiffEq.jl algorithms achieve the highest performance for non-stiff equations
while being the most generic: accepting the most Julia-based types, allow for
sophisticated event handling, etc. They are recommended for all non-stiff problems.
For stiff problems, the algorithms are currently not as high of order or as well-optimized
as the ODEInterface.jl or Sundials.jl algorithms, and thus if the problem is on
arrays of Float64, they are recommended. However, the stiff methods from OrdinaryDiffEq.jl
are able to handle a larger generality of number types (arbitrary precision, etc.)
and thus are recommended for stiff problems on non-Float64 numbers.

### Runge-Kutta Methods for Non-Stiff Equations

- `Euler`- The canonical forward Euler method. Fixed timestep only.
- `Midpoint` - The second order midpoint method. Uses embedded Euler method for
  adaptivity.
- `Heun` - The second order Heun's method. Uses embedded Euler method for
  adaptivity.
- `Ralston` - The optimized second order midpoint method. Uses embedded Euler.
  method for adaptivity.
- `RK4` - The canonical Runge-Kutta Order 4 method. Uses a defect control for
  adaptive stepping using maximum error over the whole interval.
- `BS3` - Bogacki-Shampine 3/2 method.
- `OwrenZen3` - Owren-Zennaro optimized interpolantion 3/2 method (free 3th
  order interpolant).
- `OwrenZen4` - Owren-Zennaro optimized interpolantion 4/3 method (free 4th
  order interpolant).
- `OwrenZen5` - Owren-Zennaro optimized interpolantion 5/4 method (free 5th
  order interpolant).
- `DP5` - Dormand-Prince's 5/4 Runge-Kutta method. (free 4th order interpolant).
- `Tsit5` - Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
- `TanYam7` - Tanaka-Yamashita 7 Runge-Kutta method.
- `DP8` - Hairer's 8/5/3 adaption of the Dormand-Prince Runge-Kutta method.
  (7th order interpolant).
- `TsitPap8` - Tsitouras-Papakostas 8/7 Runge-Kutta method.
- `Feagin10` - Feagin's 10th-order Runge-Kutta method.
- `Feagin12` - Feagin's 12th-order Runge-Kutta method.
- `Feagin14` - Feagin's 14th-order Runge-Kutta method.

Example usage:

```julia
alg = Tsit5()
solve(prob,alg)  
```

Additionally, the following algorithms have a lazy interpolant:

- `BS5` - Bogacki-Shampine 5/4 Runge-Kutta method. (lazy 5th order interpolant).
- `Vern6` - Verner's "Most Efficient" 6/5 Runge-Kutta method. (lazy 6th order interpolant).
- `Vern7` - Verner's "Most Efficient" 7/6 Runge-Kutta method. (lazy 7th order interpolant).
- `Vern8` - Verner's "Most Efficient" 8/7 Runge-Kutta method. (lazy 8th order interpolant)
- `Vern9` - Verner's "Most Efficient" 9/8 Runge-Kutta method. (lazy 9th order interpolant)

These methods require a few extra steps in order to compute the high order
interpolation, but these steps are only taken when the interpolation is used.
These methods when lazy assume that the parameter vector `p` will be unchanged
between the moment of the interval solving and the interpolation. If `p` is
changed in a ContinuousCallback, or in a DiscreteCallback and the continuous
solution is used after the full solution, then set `lazy=false`.

Example:

```julia
solve(prob,Vern7()) # lazy by default
solve(prob,Vern7(lazy=false))
```


### Explicit Strong-Stability Preserving Runge-Kutta Methods for Hyperbolic PDEs (Conservation Laws)

- `SSPRK22` - The two-stage, second order strong stability preserving (SSP)
  method of Shu and Osher (SSP coefficient 1, free 2nd order SSP interpolant). Fixed timestep only.
- `SSPRK33` - The three-stage, third order strong stability preserving (SSP)
  method of Shu and Osher (SSP coefficient 1, free 2nd order SSP interpolant). Fixed timestep only.
- `SSPRK53` - The five-stage, third order strong stability preserving (SSP)
  method of Ruuth (SSP coefficient 2.65, free 3rd order Hermite interpolant). Fixed timestep only.
- `SSPRK63` - The six-stage, third order strong stability preserving (SSP)
  method of Ruuth (SSP coefficient 3.518, free 3rd order Hermite interpolant). Fixed timestep only.
- `SSPRK73` - The seven-stage, third order strong stability preserving (SSP)
  method of Ruuth (SSP coefficient 4.2879, free 3rd order Hermite interpolant). Fixed timestep only.
- `SSPRK83` - The eight-stage, third order strong stability preserving (SSP)
  method of Ruuth (SSP coefficient 5.107, free 3rd order Hermite interpolant). Fixed timestep only.
- `SSPRK432` - A  3/2 adaptive strong stability preserving (SSP) method with
  five stages (SSP coefficient 2, free 2nd order SSP interpolant).
- `SSPRK932` - A  3/2 adaptive strong stability preserving (SSP) method with
  nine stages (SSP coefficient 6, free 3rd order Hermite interpolant).
- `SSPRK54` - The five-stage, fourth order strong stability preserving (SSP)
  method of Spiteri and Ruuth (SSP coefficient 1.508, 3rd order Hermite interpolant). Fixed timestep only.
- `SSPRK104` - The ten-stage, fourth order strong stability preserving method
  of Ketcheson (SSP coefficient 6, free 3rd order Hermite interpolant). Fixed timestep only.

The SSP coefficients of the methods can be queried as `ssp_coefficient(alg)`.
All explicit SSP methods take two optional arguments `SSPXY(stage_limiter!, step_limiter!)`, where
`stage_limiter!` and `step_limiter` are functions taking arguments of the form `limiter!(u, f, t)`.
Here, `u` is the new solution value (updated inplace) after an explicit Euler stage / the whole time
step , `f` the time derivative function (semidiscretisation for PDEs), and `t` the current time. These
limiters can be used to enforce physical constraints, e.g. the positivity preserving limiters of
Zhang and Shu (Zhang, Xiangxiong, and Chi-Wang Shu. "Maximum-principle-satisfying and positivity-preserving
high-order schemes for conservation laws: survey and new developments." Proceedings of the Royal Society of
London A: Mathematical, Physical and Engineering Sciences. The Royal Society, 2011.).

### Low-Storage Methods

Up to now, there are still some improvements concerning memory consumption posible, e.g. dropping the dense
output, interpolations, callbacks etc. However, some basic methods are available.

- `CarpenterKennedy2N54` - The five-stage, fourth order low-storage method of Carpenter and Kennedy
  (free 3rd order Hermite interpolant). Fixed timestep only. Designed for hyperbolic PDEs (stability properties).

### Explicit Multistep Methods

Methods using the approximation at more than one previous mesh point to determine
the approximation at the next point are called multistep methods. These methods
tend to be more efficient as the size of the system or the cost of `f` increases.

#### Adams-Bashforth Explicit Methods

These methods require a choice of `dt`.

- `AB3` - The 3-step third order multistep method. Ralston's Second Order Method
  is used to calculate starting values.
- `AB4` - The 4-step fourth order multistep method. Runge-Kutta method of order
  4 is used to calculate starting values.  
- `AB5` - The 5-step fifth order multistep method. Runge-Kutta method of order
  4 is used to calculate starting values.  
- `ABM32` - It is third order method. In `ABM32`, `AB3` works as predictor and
  Adams Moulton 2-steps method works as Corrector. Ralston's Second Order Method
  is used to calculate starting values.  
- `ABM43` - It is fourth order method. In `ABM43`, `AB4` works as predictor and
  Adams Moulton 3-steps method works as Corrector. Runge-Kutta method of order
  4 is used to calculate starting values.  
- `ABM54` - It is fifth order method. In `ABM54`, `AB5` works as predictor and
  Adams Moulton 4-steps method works as Corrector. Runge-Kutta method of order 4
  is used to calculate starting values.

#### Adaptive step size Adams explicit Methods

- `VCAB3` - The 3rd order Adams method. Bogacki-Shampine 3/2 method is used to calculate
  starting values.  
- `VCAB4` - The 4th order Adams method. Runge-Kutta 4 is used to calculate
  starting values.  
- `VCAB5` - The 5th order Adams method. Runge-Kutta 4 is used to calculate
  starting values.
- `VCABM3` - The 3rd order Adams-Moulton method. Bogacki-Shampine 3/2 method is used
  to calculate starting values.  
- `VCABM4` - The 4th order Adams-Moulton method. Runge-Kutta 4 is used to calculate
  starting values.  
- `VCABM5` - The 5th order Adams-Moulton method. Runge-Kutta 4 is used to calculate
  starting values.
- `VCABM` - An adaptive order adaptive time Adams Moulton method. It uses an
  order adaptivity algorithm is derived from Shampine's DDEABM.
- `AN5` - An adaptive 5th order fixed-leading coefficient Adams method in
  Nordsieck form.
- `JVODE_Adams` - An adaptive time adaptive order fixed-leading coefficient Adams
  method in Nordsieck form. The order adaptivity algorithm is derived from
  Sundials' `CVODE_Adams`. In development.

### Methods for Stiff Equations

#### SDIRK Methods

- `ImplicitEuler` - A 1st order implicit solver. A-B-L-stable. Adaptive
  timestepping through a divided differences estimate via memory. Strong-stability
  presurving (SSP).
- `ImplicitMidpoint` - A second order A-stable symplectic and symmetric implicit
  solver. Good for highly stiff equations which need symplectic integration.
- `Trapezoid` - A second order A-L-stable symmetric ESDIRK method. "Almost
  symplectic" without numerical dampening. Also known as Crank-Nicholson when
  applied to PDEs. Adaptive timestepping via divided differences on the memory.
  Good for highly stiff equations which are non-oscillatory.
- `TRBDF2` - A second order A-B-L-S-stable one-step ESDIRK method. Includes
  stiffness-robust error estimates for accurate adaptive timestepping, smoothed
  derivatives for highly stiff and oscillatory problems.
- `GenericImplicitEuler` - A 1st order A-B-L-stable implicit solver with adaptive
  timestepping through a divided differences estimate via memory. Strong-stability
  presurving (SSP). Uses an external nonlinear solver. Defaults to trust region
  dogleg with full Newton, making it more robust to numerical instability at
  the cost of being less efficient.
- `GenericTrapezoid` - A second order A-stable symplectic implicit solver. Also known
  as Crank-Nicholson when applied to PDEs. Adaptive timestepping via divided
  differences on the memory. Good for highly stiff equations which are non-oscillatory.
  Uses an external nonlinear solver. Defaults to trust region
  dogleg with full Newton, making it more robust to numerical instability at
  the cost of being less efficient.
- `SDIRK2` - An A-B-L stable 2nd order SDIRK method
- `Kvaerno3` - An A-L stable stiffly-accurate 3rd order ESDIRK method
- `KenCarp3` - An A-L stable stiffly-accurate 3rd order ESDIRK method with splitting
- `Cash4` - An A-L stable 4th order SDIRK method
- `Hairer4` - An A-L stable 4rd order SDIRK method
- `Hairer42` - An A-L stable 4rd order SDIRK method
- `Kvaerno4` - An A-L stable stiffly-accurate 4rd order ESDIRK method
- `KenCarp4` - An A-L stable stiffly-accurate 4rd order ESDIRK method with splitting
- `Kvaerno5` - An A-L stable stiffly-accurate 5rd order ESDIRK method
- `KenCarp5` - An A-L stable stiffly-accurate 5rd order ESDIRK method with splitting

#### Rosenbrock Methods

- `Rosenbrock23` - An Order 2/3 L-Stable Rosenbrock-W method which is good for very
  stiff equations with oscillations at low tolerances. 2nd order stiff-aware
  interpolation.
- `Rosenbrock32` - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildy
  stiff equations without oscillations at low tolerances. Note that this method
  is prone to instability in the presence of oscillations, so use with caution.
  2nd order stiff-aware interpolation.
- `ROS3P` - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high accuracy
  on discretizations of nonlinear parabolic PDEs.
- `Rodas3` - 3rd order A-stable and stiffly stable Rosenbrock method.
- `RosShamp4`- An A-stable 4th order Rosenbrock method.
- `Veldd4` - A 4th order D-stable Rosenbrock method.
- `Velds4` - A 4th order A-stable Rosenbrock method.
- `GRK4T` - An efficient 4th order Rosenbrock method.
- `GRK4A` - An A-stable 4th order Rosenbrock method. Essentially "anti-L-stable"
  but efficient.
- `Ros4LStab` - A 4th order L-stable Rosenbrock method.
- `Rodas4` - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware
  3rd order interpolant
- `Rodas42` - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware
  3rd order interpolant
- `Rodas4P` - A 4th order A-stable stiffly stable Rosenbrock method with a stiff-aware
  3rd order interpolant. 4th order on linear parabolic problems and 3rd order accurate
  on nonlinear parabolic problems (as opposed to lower if not corrected).
- `Rodas5` - A 5th order A-stable stiffly stable Rosenbrock method. Currently has
  a Hermite interpolant because its stiff-aware 3rd order interpolant is not
  yet implemented.

#### Exponential Rosenbrock Methods

- `LawsonEuler` - First order exponential Euler scheme. Fixed timestepping only.
- `NorsettEuler` - First order exponential-RK scheme. Fixed timestepping only.
  Alias: `ETD1`.
- `ETD2` - Second order Exponential Time Differencing method. Fixed timestepping only.
- `ETDRK4` - 4th order exponential-RK scheme. Fixed timestepping only.
- `HochOst4` - 4th order exponential-RK scheme with stiff order 4. Fixed
  timestepping only.

#### Exponential Propagation Iterative Runge-Kutta Methods (EPIRK)

- `Exp4` - 4th order EPIRK scheme. Fixed time stepping only.
- `EPIRK4s3A` - 4th order EPIRK scheme with stiff order 4. Fixed time stepping only.
- `EPIRK4s3B` - 4th order EPIRK scheme with stiff order 4. Fixed time stepping only.
- `EPIRK5-P1` - 5th order EPIRK scheme. Fixed time stepping only.
- `EPIRK5-P2` - 5th order EPIRK scheme. Fixed time stepping only.
- `EPIRK5s3` - 5th order "horizontal" EPIRK scheme with stiff order 5.
  Fixed time stepping only. Broken.
- `EXPRB53s3`- 5th order EPIRK scheme with stiff order 5. Fixed time stepping only.

#### Multistep Methods

- `QNDF1` - An adaptive order 1 quasi-constant timestep L-stable numerical
  differentiation function (NDF) method. Optional parameter `kappa` defaults
  to Shampine's accuracy-optimal `-0.1850`.
- `QBDF1` - An adaptive order 1 quasi-constant timestep L-stable BDF method.
  This is equivalent to implicit Euler but using the BDF error estimator with
  quasi-constant stepping.
- `ABDF2` - An adaptive order 2 L-stable fixed leading coefficient multistep
  BDF method.
- `QNDF` - An adaptive order quasi-fixed time step NDF method. Utilizes
  Shampine's accuracy-optimal `kappa` values as defaults (has a keyword argument
  for a tuple of `kappa` coefficients). In development.
- `JVODE_BDF` - An adaptive time adaptive order fixed-leading coefficient BDF
  method in Nordsieck form. The order adaptivity algorithm is derived from
  Sundials' `CVODE_BDF`. In development.

#### Implicit Strong-Stability Preserving Runge-Kutta Methods for Hyperbolic PDEs (Conservation Laws)

- `SSPSDIRK2` - A second order A-L stable symplectic SDIRK method with the strong
  stability preserving (SSP) property (SSP coefficient 2). Fixed timestep only.

#### Extra Options

All of the Rosenbrock and SDIRK methods allow for specification of `linsolve`:
the linear solver which is used. For more information on specifying the linear
solver, see
[the manual page on solver specification](../features/linear_nonlinear.html).

The following methods allow for specification of `nlsolve`: the nonlinear
solver which is used:

- `GenericImplicitEuler`
- `GenericTrapezoid`

Note that performance overload information (Jacobians etc.) are not used in this
mode. This can control autodifferentiation of the Jacobian as well.
For more information on specifying the nonlinear solver, see
[the manual page on solver specification](../features/linear_nonlinear.html).

Additionally, the Rosenbrock and SDIRK methods have differentiation
controls. In each of these, `autodiff` can be set to turn on/off
autodifferentiation, and `chunk_size` can be used to set the chunksize of the Dual
 numbers (see the
[documentation for ForwardDiff.jl for details](http://www.juliadiff.org/ForwardDiff.jl/advanced_usage.html#configuring-chunk-size)).
In addition, the Rosenbrock and SDIRK methods can set `diff_type`, which is the
type of numerical differentiation that is used (when autodifferentiation is
disabled). The choices are `Val{:central}`, `Val{:forward}` or `Val{:complex}`.

Examples:

```julia
sol = solve(prob,Rosenbrock23()) # Standard, uses autodiff
sol = solve(prob,Rosenbrock23(chunk_size=10)) # Autodiff with chunksize of 10
sol = solve(prob,Rosenbrock23(autodiff=false)) # Numerical differentiation with central differencing
sol = solve(prob,Rosenbrock23(autodiff=false,diff_type=Val{:forward})) # Numerical differentiation with forward differencing
```

### Tableau Method

Additionally, there is the tableau method:

  - `ExplicitRK` - A general Runge-Kutta solver which takes in a tableau. Can be adaptive. Tableaus
  are specified via the keyword argument `tab=tableau`. The default tableau is
  for Dormand-Prince 4/5. Other supplied tableaus can be found in the Supplied Tableaus section.

Example usage:

```julia
alg = ExplicitRK(tableau=constructDormandPrince())
solve(prob,alg)
```

### CompositeAlgorithm

One unique feature of OrdinaryDiffEq.jl is the `CompositeAlgorithm`, which allows
you to, with very minimal overhead, design a multimethod which switches between
chosen algorithms as needed. The syntax is `CompositeAlgorithm(algtup,choice_function)`
where `algtup` is a tuple of OrdinaryDiffEq.jl algorithms, and `choice_function`
is a function which declares which method to use in the following step. For example,
we can design a multimethod which uses `Tsit5()` but switches to `Vern7()` whenever
`dt` is too small:

```julia
choice_function(integrator) = (Int(integrator.dt<0.001) + 1)
alg_switch = CompositeAlgorithm((Tsit5(),Vern7()),choice_function)
```

The `choice_function` takes in an `integrator` and thus all of the features
available in the [Integrator Interface](@ref) can be used in the choice
function.

A helper algorithm was created for building 2-method automatic switching for
stiffness detection algorithms. This is the `AutoSwitch` algorithm with the
following options:

```julia
AutoSwitch(nonstiffalg::nAlg, stiffalg::sAlg;
           maxstiffstep=10, maxnonstiffstep=3,
           nonstifftol::T=9//10, stifftol::T=9//10,
           dtfac=2.0, stiffalgfirst=false)
```

The `nonstiffalg` must have an appropriate stiffness estimate built into the
method. The `stiffalg` can receive its estimate from the Jacobian calculation.
`maxstiffstep` is the number of stiffness detects before switching to the stiff
algorithm and `maxnonstiffstep` is vice versa. `nonstifftol` and `stifftol` are
the tolerances associated with the stiffness comparison against the stability
region. `dtfac` is the factor that `dt` is changed when switching: multiplied
when going from non-stiff to stiff and divided when going stiff to non-stiff.
`stiffalgfirst` denotes whether the first step should use the stiff algorithm.

#### Pre-Built Stiffness Detecting and Auto-Switching Algorithms

These methods require a `Autoalg(stiffalg)` to be chosen as the method to switch
to when the ODE is stiff. It can be any of the OrdinaryDiffEq.jl one-step stiff
methods and has all of the arguments of the `AutoSwitch` algorithm.

- `AutoTsit5` - `Tsit5` with automated switching.
- `AutoDP5` - `DP5` with automated switching.
- `AutoVern6` - `Vern6` with automated switching.
- `AutoVern7` - `Vern7` with automated switching.
- `AutoVern8` - `Vern8` with automated switching.
- `AutoVern9` - `Vern9` with automated switching.

Example:

```julia
tsidas_alg = AutoTsit5(Rodas5())
sol = solve(prob,tsidas_alg)

tsidas_alg = AutoTsit5(Rodas5(),nonstifftol = 11/10)
```

Is the `Tsit5` method with automatic switching to `Rodas5`.

## Sundials.jl

The Sundials suite is built around multistep methods. These methods are more efficient
than other methods when the cost of the function calculations is really high, but
for less costly functions the cost of nurturing the timestep overweighs the benefits.
However, the BDF method is a classic method for stiff equations and "generally works".

  - `CVODE_BDF` - CVode Backward Differentiation Formula (BDF) solver.
  - `CVODE_Adams` - CVode Adams-Moulton solver.
  - `ARKODE` - Explicit and ESDIRK Runge-Kutta methods of orders 2-8 depending
    on choice of options.

The Sundials algorithms all come with a 3rd order Hermite polynomial interpolation.
Note that the constructors for the Sundials algorithms take two main arguments:

  - `method` - This is the method for solving the implicit equation. For BDF this
    defaults to `:Newton` while for Adams this defaults to `:Functional`. These
    choices match the recommended pairing in the Sundials.jl manual. However,
    note that using the `:Newton` method may take less iterations but requires
    more memory than the `:Function` iteration approach.
  - `linearsolver` - This is the linear solver which is used in the `:Newton` method.

  The choices for the linear solver are:

  - `:Dense` - A dense linear solver.
  - `:Band` - A solver specialized for banded Jacobians. If used, you must set the
    position of the upper and lower non-zero diagonals via `jac_upper` and
    `jac_lower`.
  - `:Diagonal` - This method is specialized for diagonal Jacobians.
  - `:GMRES` - A GMRES method. Recommended first choice Krylov method
  - `:BCG` - A Biconjugate gradient method.
  - `:PCG` - A preconditioned conjugate gradient method. Only for symmetric linear systems.
  - `:TFQMR` - A TFQMR method.
  - `:KLU` - A sparse factorization method. Requires that the user specifies a
    Jacobian. The Jacobian must be set as a sparse matrix in the `ODEProblem`
    type.

Example:

```julia
CVODE_BDF() # BDF method using Newton + Dense solver
CVODE_BDF(method=:Functional) # BDF method using Functional iterations
CVODE_BDF(linear_solver=:Band,jac_upper=3,jac_lower=3) # Banded solver with nonzero diagonals 3 up and 3 down
CVODE_BDF(linear_solver=:BCG) # Biconjugate gradient method                                   
```

The main options for `ARKODE` are the choice between explicit and implicit and
the method order, given via:

```julia
ARKODE(Sundials.Explicit()) # Solve with explicit tableau of default order 4
ARKODE(Sundials.Implicit(),order = 3) # Solve with explicit tableau of order 3
```

The order choices for explicit are 2 through 8 and for implicit 3 through 5.
Specific methods can also be set through the `etable` and `itable` options
for explicit and implicit tableaus respectively. The available tableaus are:

`etable`:

- `HEUN_EULER_2_1_2`: 2nd order Heun's method
- `BOGACKI_SHAMPINE_4_2_3`:
- `ARK324L2SA_ERK_4_2_3`: explicit portion of Kennedy and Carpenter's 3rd
  order method
- `ZONNEVELD_5_3_4`: 4th order explicit method
- `ARK436L2SA_ERK_6_3_4`: explicit portion of Kennedy and Carpenter's 4th
  order method
- `SAYFY_ABURUB_6_3_4`: 4th order explicit method
- `CASH_KARP_6_4_5`: 5th order explicit method
- `FEHLBERG_6_4_5`: Fehlberg's classic 5th order method
- `DORMAND_PRINCE_7_4_5`: the classic 5th order Dormand-Prince method
- `ARK548L2SA_ERK_8_4_5`: explicit portion of Kennedy and Carpenter's 5th
  order method
- `VERNER_8_5_6`: Verner's classic 5th order method
- `FEHLBERG_13_7_8`: Fehlberg's 8th order method

`itable`:

- `SDIRK_2_1_2`: An A-B-stable 2nd order SDIRK method
- `BILLINGTON_3_3_2`: A second order method with a 3rd order error predictor
  of less stability
- `TRBDF2_3_3_2`: The classic TR-BDF2 method
- `KVAERNO_4_2_3`: an L-stable 3rd order ESDIRK method
- `ARK324L2SA_DIRK_4_2_3`: implicit portion of Kennedy and Carpenter's 3th
  order method
- `CASH_5_2_4`: Cash's 4th order L-stable SDIRK method
- `CASH_5_3_4`: Cash's 2nd 4th order L-stable SDIRK method
- `SDIRK_5_3_4`: Hairer's 4th order SDIRK method
- `KVAERNO_5_3_4`: Kvaerno's 4th order ESDIRK method
- `ARK436L2SA_DIRK_6_3_4`: implicit portion of Kennedy and Carpenter's 4th
  order method
- `KVAERNO_7_4_5`: Kvaerno's 5th order ESDIRK method
- `ARK548L2SA_DIRK_8_4_5`: implicit portion of Kennedy and Carpenter's 5th
  order method

These can be set for example via:

```julia
ARKODE(Sundials.Explicit(),etable = Sundials.DORMAND_PRINCE_7_4_5)
ARKODE(Sundials.Implicit(),itable = Sundials.KVAERNO_4_2_3)
```

All of the additional options are available. The full constructor is:

```julia
CVODE_BDF(;method=:Newton,linear_solver=:Dense,
          jac_upper=0,jac_lower=0,
          stored_upper = jac_upper + jac_lower,
          non_zero=0,krylov_dim=0,
          stability_limit_detect=false,
          max_hnil_warns = 10,
          max_order = 5,
          max_error_test_failures = 7,
          max_nonlinear_iters = 3,
          max_convergence_failures = 10)

CVODE_Adams(;method=:Functional,linear_solver=:None,
            jac_upper=0,jac_lower=0,
            stored_upper = jac_upper + jac_lower,
            krylov_dim=0,
            stability_limit_detect=false,
            max_hnil_warns = 10,
            max_order = 12,
            max_error_test_failures = 7,
            max_nonlinear_iters = 3,
            max_convergence_failures = 10)

ARKODE(stiffness=Sundials.Implicit();
      method=:Newton,linear_solver=:Dense,
      jac_upper=0,jac_lower=0,stored_upper = jac_upper+jac_lower,
      non_zero=0,krylov_dim=0,
      max_hnil_warns = 10,
      max_error_test_failures = 7,
      max_nonlinear_iters = 3,
      max_convergence_failures = 10,
      predictor_method = 0,
      nonlinear_convergence_coefficient = 0.1,
      dense_order = 3,
      order = 4,
      set_optimal_params = false,
      crdown = 0.3,
      dgmax = 0.2,
      rdiv = 2.3,
      msbp = 20,
      adaptivity_method = 0
      )
```

See [the CVODE manual](https://computation.llnl.gov/sites/default/files/public/cv_guide.pdf)
and the [ARKODE manual](https://computation.llnl.gov/sites/default/files/public/ark_guide.pdf)
for details on the additional options.

## ODEInterface.jl

The ODEInterface algorithms are the classic Fortran algorithms. While the
non-stiff algorithms are superseded by the more featured and higher performance
Julia implementations from OrdinaryDiffEq.jl, the stiff solvers such as `radau`
are some of the most efficient methods available (but are restricted for use on
arrays of Float64).

Note that this setup is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use ODEInterfaceDiffEq.jl:

```julia
Pkg.add("ODEInterfaceDiffEq")
using ODEInterfaceDiffEq
```

  - `dopri5` - Hairer's classic implementation of the Dormand-Prince 4/5 method.
  - `dop853` - Explicit Runge-Kutta 8(5,3) by Dormand-Prince.
  - `odex` - GBS extrapolation-algorithm based on the midpoint rule.
  - `seulex` - Extrapolation-algorithm based on the linear implicit Euler method.
  - `radau` - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.
  - `radau5` - Implicit Runge-Kutta method (Radau IIA) of order 5.
  - `rodas` - Rosenbrock 4(3) method.
  - `ddeabm` - Adams-Bashforth-Moulton Predictor-Corrector method (order between 1 and 12)
  - `ddebdf` - Backward Differentiation Formula (orders between 1 and 5)

Note that while the output only has a linear interpolation, a higher order
interpolation is used for intermediate dense output for `saveat` and for
event handling.

## LSODA.jl

This setup provides a wrapper to the algorithm LSODA, a well-known method which uses switching
to solve both stiff and non-stiff equations.

  - `lsoda` - The LSODA wrapper algorithm.

Note that this setup is not automatically included with DifferentialEquaitons.jl.
To use the following algorithms, you must install and use LSODA.jl:

```julia
Pkg.add("LSODA")
using LSODA
```

## ODE.jl

  - `ode23` - Bogacki-Shampine's order 2/3 Runge-Kutta  method
  - `ode45` - A Dormand-Prince order 4/5 Runge-Kutta method
  - `ode23s` - A modified Rosenbrock order 2/3 method due to Shampine
  - `ode78` - A Fehlburg order 7/8 Runge-Kutta method
  - `ode4` - The classic Runge-Kutta order 4 method
  - `ode4ms` - A fixed-step, fixed order Adams-Bashforth-Moulton method†
  - `ode4s` - A 4th order Rosenbrock method due to Shampine


†: Does not step to the interval endpoint. This can cause issues with discontinuity
detection, and [discrete variables need to be updated appropriately](../features/diffeq_arrays.html).

## MATLABDiffEq.jl

These algorithms require that the problem was defined using a `ParameterizedFunction`
via the `@ode_def` macro. Note that this setup is not automatically included
with DifferentialEquaitons.jl. To use the following algorithms, you must install
and use MATLABDiffEq.jl:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/MATLABDiffEq.jl")
using MATLABDiffEq
```

This requires a licensed MATLAB installation. The available methods are:

  - `ode23`
  - `ode45`
  - `ode113`
  - `ode23s`
  - `ode23t`
  - `ode23tb`
  - `ode15s`
  - `ode15i`

For more information on these algorithms, see
[the MATLAB documentation](https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html).

## GeometricIntegrators.jl

GeometricIntegrators.jl is a set of fixed timestep algorithms written in Julia.
Note that this setup is not automatically included with DifferentialEquaitons.jl.
To use the following algorithms, you must install and use
GeometricIntegratorsDiffEq.jl:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/GeometricIntegratorsDiffEq.jl")
using GeometricIntegratorsDiffEq
```

- `GIEuler` - 1st order Euler method
- `GIMidpoint` - 2nd order explicit midpoint method
- `GIHeun` - 2nd order Heun's method
- `GIKutta` - 3rd order Kutta's method
- `GIERK4` - standard 4th order Runge-Kutta
- `GIERK438` - 4th order Runge-Kutta, 3/8's rule
- `GIImplicitEuler` - 1st order implicit Euler method
- `GIImplicitMidpoint` - 2nd order implicit midpoint method
- `GIRadIIA2` - 2-stage order 3 Radau-IIA
- `GIRadIIA3` - 3-stage order 5 Radau-IIA
- `GISRK3` - 3-stage order 4 symmetric Runge-Kutta method
- `GIGLRK(s)` - Gauss-Legendre Runge-Kutta method of order 2s

Note that all of these methods require the user supplies `dt`.

## BridgeDiffEq.jl

Bridge.jl is a set of fixed timestep algorithms written in Julia. These methods
are made and optimized for out-of-place functions on immutable (static vector)
types. Note that this setup is not automatically included with
DifferentialEquaitons.jl. To use the following algorithms, you must install and
use BridgeDiffEq.jl:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/BridgeDiffEq.jl")
using BridgeDiffEq
```

- `BridgeR3` - 3rd order Ralston method
- `BridgeBS3` - 3rd order Bogacki-Shampine method

## TaylorIntegration.jl

TaylorIntegration.jl is a pure-Julia implementation of an adaptive order Taylor
series method for high accuracy integration of ODEs. These methods are optimized
when the absolute tolerance is required to be very low.
Note that this setup is not automatically included with DifferentialEquaitons.jl.
To use the following algorithms, you must install and
use TaylorIntegration.jl:

```julia
Pkg.add("TaylorIntegration")
using TaylorIntegration
```

- `TaylorMethod(order)` - Taylor integration method with maximal `order` (required)

## List of Supplied Tableaus

A large variety of tableaus have been supplied by default via DiffEqDevTools.jl.
The list of tableaus can be found in [the developer docs](https://juliadiffeq.github.io/DiffEqDevDocs.jl/latest/internals/tableaus.html).
For the most useful and common algorithms, a hand-optimized version is supplied
in OrdinaryDiffEq.jl which is recommended for general uses (i.e. use `DP5`
instead of `ExplicitRK` with `tableau=constructDormandPrince()`). However, these
serve as a good method for comparing between tableaus and understanding the
pros/cons of the methods. Implemented are every published tableau (that I know
exists). Note that user-defined tableaus also are accepted. To see how to define
a tableau, checkout the [premade tableau source code](https://github.com/JuliaDiffEq/DiffEqDevTools.jl/blob/master/src/ode_tableaus.jl).
Tableau docstrings should have appropriate citations (if not, file an issue).

Plot recipes are provided which will plot the stability region for a given tableau.
