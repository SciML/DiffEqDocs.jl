# ODE Solvers

`solve(prob::ODEProblem,alg;kwargs)`

Solves the ODE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

It is suggested that you try choosing an algorithm using the `alg_hints`
keyword argument. However, in some cases you may want something specific,
or you may just be curious. This guide is to help you choose the right algorithm.

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
very large (`>10,000` ODEs?), the function calculation is very expensive,
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
systems `TRBDF2` will be more efficient.

At medium tolerances (`>1e-8`?) it is recommended you use `Rodas5`,
`Rodas4P` (the former is more efficient but the later is more reliable),
`Kvaerno5`, or `KenCarp4`. As native DifferentialEquations.jl solvers,
many Julia numeric types (such as BigFloats,
[ArbFloats](https://github.com/JuliaArbTypes/ArbFloats.jl), or
[DecFP](https://github.com/stevengj/DecFP.jl)) will work. When the equation is
defined via the `@ode_def` macro, these will be the most efficient.

For faster solving at low tolerances (`<1e-9`) but when `Vector{Float64}` is used,
use `radau`.

For asymtopically large systems of ODEs (`N>10000`?)
where `f` is very costly and the complex eigenvalues are minimal (low oscillations),
in that case `CVODE_BDF` will be the most efficient but requires `Vector{Float64}`.
`CVODE_BDF` will also do surprisingly well if the solution is smooth.

#### Special Properties of Stiff Integrators

`Trapezoid` is a symplectic integrator. `ImplicitEuler` is an extension
to the common algorithm with adaptive timestepping and efficient quasi-Newton
Jacobian reusage which is strong-stability presurving (SSP) for hyperbolic PDEs.

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
- `RK4` - The canonical Runge-Kutta Order 4 method. Uses a defect control to
  adaptive step using maximum error over the whole interval.
- `BS3` - Bogacki-Shampine 3/2 method.
- `OwrenZen3` - Owren-Zennaro optimized interpolantion 3/2 method (free 3th order interpolant).
- `OwrenZen4` - Owren-Zennaro optimized interpolantion 4/3 method (free 4th order interpolant).
- `OwrenZen5` - Owren-Zennaro optimized interpolantion 5/4 method (free 5th order interpolant).
- `DP5` - Dormand-Prince's 5/4 Runge-Kutta method. (free 4th order interpolant).
- `Tsit5` - Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
- `BS5` - Bogacki-Shampine 5/4 Runge-Kutta method. (5th order interpolant).
- `Vern6` - Verner's "Most Efficient" 6/5 Runge-Kutta method. (6th order interpolant).
- `Vern7` - Verner's "Most Efficient" 7/6 Runge-Kutta method. (7th order interpolant).
- `TanYam7` - Tanaka-Yamashita 7 Runge-Kutta method.
- `DP8` - Hairer's 8/5/3 adaption of the Dormand-Prince 8
  method Runge-Kutta method. (7th order interpolant).
- `TsitPap8` - Tsitouras-Papakostas 8/7 Runge-Kutta method.
- `Vern8` - Verner's "Most Efficient" 8/7 Runge-Kutta method. (8th order interpolant)
- `Vern9` - Verner's "Most Efficient" 9/8 Runge-Kutta method. (9th order interpolant)
- `Feagin10` - Feagin's 10th-order Runge-Kutta method.
- `Feagin12` - Feagin's 12th-order Runge-Kutta method.
- `Feagin14` - Feagin's 14th-order Runge-Kutta method.

Example usage:

```julia
alg = Tsit5()
solve(prob,alg)  
```

### Explicit Strong-Stability Presurving Runge-Kutta Methods for Hyperbolic PDEs (Conservation Laws)

- `SSPRK22` - The two-stage, second order strong stability preserving (SSP)
  method of Shu and Osher. (free 2nd order SSP interpolant). Fixed timestep only.
- `SSPRK33` - The three-stage, third order strong stability preserving (SSP)
  method of Shu and Osher. (free 2nd order SSP interpolant). Fixed timestep only.
- `SSPRK432` - A  3/2 adaptive strong stability preserving (SSP) method with
  five stages. (free 2nd order SSP interpolant)
- `SSPRK104` - The ten-stage, fourth order strong stability preserving method
  of Ketcheson. (free 3rd order Hermite interpolant). Fixed timestep only.

### Methods for Stiff Equations

#### SDIRK Methods

- `ImplicitEuler` - A 1st order implicit solver. A-B-L-stable. Adaptive
  timestepping through a divided differences estimate via memory. Strong-stability
  presurving (SSP).
- `Trapezoid` - A second order A-stable symplectic implicit solver. Also known
  as Crank-Nicholson when applied to PDEs. Adaptive timestepping via divided
  differences on the memory. Good for highly stiff equations which are non-oscillatory.
- `TRBDF2` - A second order A-B-L-S-stable one-step method. Includes
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
- `Kvaerno3` - An A-L stable stiffly-accurate 3rd order SDIRK method
- `KenCarp3` - An A-L stable stiffly-accurate 3rd order SDIRK method with splitting
- `Cash4` - An A-L stable 4th order SDIRK method
- `Hairer4` - An A-L stable 4rd order SDIRK method
- `Hairer42` - An A-L stable 4rd order SDIRK method
- `Kvaerno4` - An A-L stable stiffly-accurate 4rd order SDIRK method
- `KenCarp4` - An A-L stable stiffly-accurate 4rd order SDIRK method with splitting
- `Kvaerno5` - An A-L stable stiffly-accurate 5rd order SDIRK method
- `KenCarp5` - An A-L stable stiffly-accurate 5rd order SDIRK method with splitting

#### Rosenbrock Methods

- `Rosenbrock23` - An Order 2/3 L-Stable Rosenbrock-W method which is good for very
  stiff equations with oscillations at low tolerances.
- `Rosenbrock32` - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildy
  stiff equations without oscillations at low tolerances. Note that this method
  is prone to instability in the presence of oscillations, so use with caution.
- `ROS3P` - 3rd order A-stable and stiffly stable (Index-1 DAE compatible)
  Rosenbrock method. Keeps high accuracy on discretizations of nonlinear parabolic
  PDEs.
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
- `Rodas5` - A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware
  3rd order interpolant.

### Implicit Strong-Stability Presurving Runge-Kutta Methods for Hyperbolic PDEs (Conservation Laws)

- `SSPSDIRK2` - A second order A-L stable symplectic SDIRK method with the strong
  stability presurving (SSP) property.

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
disabled). The choices are `:central` or `:forward`.

Examples:

```julia
sol = solve(prob,Rosenbrock23()) # Standard, uses autodiff
sol = solve(prob,Rosenbrock23(chunk_size=10)) # Autodiff with chunksize of 10
sol = solve(prob,Rosenbrock23(autodiff=false)) # Numerical differentiation with central differencing
sol = solve(prob,Rosenbrock23(autodiff=false,diff_type=:forward)) # Numerical differentiation with forward differencing
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
chosen algorithms as needed. The syntax is `CompositeAlgorthm(algtup,choice_function)`
where `algtup` is a tuple of OrdinaryDiffEq.jl algorithms, and `choice_function`
is a function which declares which method to use in the following step. For example,
we can design a multimethod which uses `Tsit5()` but switches to `Vern7()` whenever
`dt` is too small:

```julia
choice_function(integrator) = (Int(integrator.dt<0.001) + 1)
alg_switch = CompositeAlgorithm((Tsit5(),Vern7()),choice_function)
```

The `choice_function` takes in an `integrator` and thus all of the features
available in the [Integrator Interface](@ref)
can be used in the choice function.

## Sundials.jl

The Sundials suite is built around multistep methods. These methods are more efficient
than other methods when the cost of the function calculations is really high, but
for less costly functions the cost of nurturing the timestep overweighs the benefits.
However, the BDF method is a classic method for stiff equations and "generally works".

  - `CVODE_BDF` - CVode Backward Differentiation Formula (BDF) solver.
  - `CVODE_Adams` - CVode Adams-Moulton solver.

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
  - `:BCG` - A Biconjugate gradient method.
  - `:TFQMR` - A TFQMR method.

Example:

```julia
CVODE_BDF() # BDF method using Newton + Dense solver
CVODE_BDF(method=:Functional) # BDF method using Functional iterations
CVODE_BDF(linear_solver=:Band,jac_upper=3,jac_lower=3) # Banded solver with nonzero diagonals 3 up and 3 down
CVODE_BDF(linear_solver=:BCG) # Biconjugate gradient method                                   
```

All of the additional options are available. The full constructor is:

```julia
CVODE_BDF(;method=:Newton,linear_solver=:Dense,
          jac_upper=0,jac_lower=0,non_zero=0,krylov_dim=0,
          stability_limit_detect=false,
          max_hnil_warns = 10,
          max_order = 5,
          max_error_test_failures = 7,
          max_nonlinear_iters = 3,
          max_convergence_failures = 10)

CVODE_Adams(;method=:Functional,linear_solver=:None,
            jac_upper=0,jac_lower=0,krylov_dim=0,
            stability_limit_detect=false,
            max_hnil_warns = 10,
            max_order = 12,
            max_error_test_failures = 7,
            max_nonlinear_iters = 3,
            max_convergence_failures = 10)
```

See [the Sundials manual](https://computation.llnl.gov/sites/default/files/public/cv_guide.pdf)
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

## List of Supplied Tableaus

A large variety of tableaus have been supplied by default via DiffEqDevTools.jl.
The list of tableaus can be found in [the developer docs](https://juliadiffeq.github.io/DiffEqDevDocs.jl/latest/internals/tableaus.html).
For the most useful and common algorithms, a hand-optimized version is supplied
in OrdinaryDiffEq.jl which is recommended for general uses (i.e. use `DP5` instead of `ExplicitRK`
with `tableau=constructDormandPrince()`). However, these serve as a good method
for comparing between tableaus and understanding the pros/cons of the methods.
Implemented are every published tableau (that I know exists). Note that user-defined
tableaus also are accepted. To see how to define a tableau, checkout the [premade tableau source code](https://github.com/JuliaDiffEq/DiffEqDevTools.jl/blob/master/src/ode_tableaus.jl).
Tableau docstrings should have appropriate citations (if not, file an issue).

Plot recipes are provided which will plot the stability region for a given tableau.
