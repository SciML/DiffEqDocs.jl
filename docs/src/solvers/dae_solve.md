# Mass Matrix and Fully Implicit DAE Solvers

## Recommended Methods

For medium to low accuracy small numbers of DAEs in constant mass matrix form,
the  `Rosenbrock23` and `Rodas4` methods are good choices which will get good
efficiency if the mass matrix is constant. `Rosenbrock23` is better for low
accuracy (error tolerance `<1e-4`) and `Rodas4` is better for high accuracy.
Another choice at high accuracy is `Rodas5P` and `RadauIIA5`.

Non-constant mass matrices are not directly supported: users are advised to
transform their problem through substitution to a DAE with constant mass
matrices.

If the problem cannot be defined in mass matrix form, the recommended method for
performance is `IDA` from the Sundials.jl package if you are solving problems with
`Float64`. If Julia types are required, currently `DFBDF` is the best method
but still needs more optimizations.

## [Full List of Methods](@id dae_solve_full)

#### Initialization Schemes

For all OrdinaryDiffEq.jl methods, an initialization scheme can be set with a
common keyword argument `initializealg`. The choices are:

  - `CheckInit`: Check that the provided initial conditions satisfy the equation;
    if not, error. This avoids the need for nonlinear solution, as well as avoiding
    changing the provided initial conditions.
  - `BrownFullBasicInit`: For Index-1 DAEs implicit DAEs and semi-explicit
    DAEs in mass matrix form. Keeps the differential variables constant. Requires
    `du0` when used on a `DAEProblem`.
  - `ShampineCollocationInit`: For Index-1 DAEs implicit DAEs and semi-explicit
    DAEs in mass matrix form. Changes both the differential and algebraic variables.
  - `NoInit`: Explicitly opts-out of DAE initialization.

### OrdinaryDiffEq.jl (Implicit ODE)

These methods from OrdinaryDiffEq are for `DAEProblem` specifications.

  - `DImplicitEuler` - 1st order A-L and stiffly stable adaptive implicit Euler
  - `DABDF2` - 2nd order A-L stable adaptive BDF method.
  - `DFBDF` - A fixed-leading coefficient adaptive-order adaptive-time BDF method,
    similar to `ode15i` or `IDA` in divided differences form.

### OrdinaryDiffEq.jl (Mass Matrix)

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

!!! note
    
    The standard Hermite interpolation used for ODE methods in OrdinaryDiffEq.jl
    falls back to a linear interpolation on the differential variables. If the
    mass matrix is non-diagonal, the Hermite interpolation does not have a fallback
    and will error.

#### Rosenbrock Methods

  - `ROS3P` - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high
    accuracy on discretizations of nonlinear parabolic PDEs.
  - `Rodas3` - 3rd order A-stable and stiffly stable Rosenbrock method.
  - `RosShamp4`- An A-stable 4th order Rosenbrock method.
  - `Veldd4` - A 4th order D-stable Rosenbrock method.
  - `Velds4` - A 4th order A-stable Rosenbrock method.
  - `GRK4T` - An efficient 4th order Rosenbrock method.
  - `GRK4A` - An A-stable 4th order Rosenbrock method. Essentially "anti-L-stable"
    but efficient.
  - `Ros4LStab` - A 4th order L-stable Rosenbrock method.
  - `Rodas4` - A 4th order A-stable stiffly stable Rosenbrock method with a
    stiff-aware 3rd order interpolant
  - `Rodas42` - A 4th order A-stable stiffly stable Rosenbrock method with a
    stiff-aware 3rd order interpolant
  - `Rodas4P` - A 4th order A-stable stiffly stable Rosenbrock method with a
    stiff-aware 3rd order interpolant. 4th order on linear parabolic problems
    and 3rd order accurate on nonlinear parabolic problems (as opposed to lower
    if not corrected).
  - `Rodas4P2` - A 4th order L-stable stiffly stable Rosenbrock method with a stiff-aware
    3rd order interpolant. 4th order on linear parabolic problems and 3rd order accurate
    on nonlinear parabolic problems. It is an improvement of Roadas4P and in case of
    inexact Jacobians a second order W method.
  - `Rodas5` - A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware
    4th order interpolant.
  - `Rodas5P` - A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware
    4th order interpolant. Has improved stability in the adaptive time stepping embedding.

#### Rosenbrock-W Methods

  - `Rosenbrock23` - An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.
  - `Rosenbrock32` - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildly stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.
  - `RosenbrockW6S4OS` - A 4th order L-stable Rosenbrock-W method (fixed step only).
  - `ROS34PW1a` - A 4th order L-stable Rosenbrock-W method.
  - `ROS34PW1b` - A 4th order L-stable Rosenbrock-W method.
  - `ROS34PW2` - A 4th order stiffy accurate Rosenbrock-W method for PDAEs.
  - `ROS34PW3` - A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method.

!!! note
    
    `Rosenbrock23` and `Rosenbrock32` have a stiff-aware interpolation but this interpolation is not safe for the algebraic variables.
    Thus use the interpolation (and therefore `saveat`) with caution if the default Hermite interpolation is used.

#### FIRK Methods

  - `RadauIIA5` - An A-B-L stable fully implicit Runge-Kutta method with internal
    tableau complex basis transform for efficiency.

#### SDIRK Methods

  - `ImplicitEuler` - Stage order 1. A-B-L-stable. Adaptive
    timestepping through a divided differences estimate via memory. Strong-stability
    preserving (SSP).
  - `ImplicitMidpoint` - Stage order 1. Symplectic. Good for when symplectic
    integration is required.
  - `Trapezoid` - A second order A-stable symmetric ESDIRK method. "Almost
    symplectic" without numerical dampening. Also known as Crank-Nicolson when
    applied to PDEs. Adaptive timestepping via divided differences on the memory.
    Good for highly stiff equations which are non-oscillatory.

#### Multistep Methods

Quasi-constant stepping is the time stepping strategy which matches the classic
GEAR, LSODE,  and `ode15s` integrators. The variable-coefficient methods match
the ideas of the classic EPISODE integrator and early VODE designs. The Fixed
Leading Coefficient (FLC) methods match the behavior of the classic VODE and
Sundials CVODE integrator.

  - `QNDF1` - An adaptive order 1 quasi-constant timestep L-stable numerical
    differentiation function (NDF) method. Optional parameter `kappa` defaults
    to Shampine's accuracy-optimal `-0.1850`.
  - `QBDF1` - An adaptive order 1 L-stable BDF method. This is equivalent to
    implicit Euler but using the BDF error estimator.
  - `ABDF2` - An adaptive order 2 L-stable fixed leading coefficient multistep
    BDF method.
  - `QNDF2` - An adaptive order 2 quasi-constant timestep L-stable numerical
    differentiation function (NDF) method.
  - `QBDF2` - An adaptive order 2 L-stable BDF method using quasi-constant timesteps.
  - `QNDF` - An adaptive order quasi-constant timestep NDF method. Utilizes
    Shampine's accuracy-optimal `kappa` values as defaults (has a keyword argument
    for a tuple of `kappa` coefficients).
  - `QBDF` - An adaptive order quasi-constant timestep BDF method.
  - `FBDF` - A fixed-leading coefficient adaptive-order adaptive-time BDF method,
    similar to `ode15i` or `CVODE_BDF` in divided differences form.

### [Sundials.jl](@id dae_solve_sundials)

Note that this setup is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use Sundials.jl:

```julia
using Pkg
Pkg.add("Sundials")
using Sundials
```

  - `IDA`: A fixed-leading coefficient fully implicit BDF method. Efficient for large systems.

For more details on controlling the Sundials.jl solvers, see the
[Sundials detailed solver API page](@ref sundials)

### DASKR.jl

DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
using Pkg
Pkg.add("DASKR")
using DASKR
```

  - `daskr` - This is a wrapper for the well-known DASKR algorithm.

For more details on controlling the DASKR.jl solvers, see the
[DASKR detailed solver API page](@ref daskr)

### DASSL.jl

  - `dassl` - A native Julia implementation of the DASSL algorithm.

### ODEInterfaceDiffEq.jl

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

  - `seulex` - Extrapolation-algorithm based on the linear implicit Euler method.
  - `radau` - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.
  - `radau5` - Implicit Runge-Kutta method (Radau IIA) of order 5.
  - `rodas` - Rosenbrock 4(3) method.
