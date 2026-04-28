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

  - `CheckInit` *(from `SciMLBase`)*: Check that the provided initial conditions satisfy the equation;
    if not, error. This avoids the need for nonlinear solution, as well as avoiding
    changing the provided initial conditions.
  - `BrownFullBasicInit` *(from `OrdinaryDiffEqNonlinearSolve`)*: For Index-1 DAEs implicit DAEs and semi-explicit
    DAEs in mass matrix form. Keeps the differential variables constant. Requires
    `du0` when used on a `DAEProblem`.
  - `ShampineCollocationInit` *(from `OrdinaryDiffEqNonlinearSolve`)*: For Index-1 DAEs implicit DAEs and semi-explicit
    DAEs in mass matrix form. Changes both the differential and algebraic variables.
  - `NoInit` *(from `SciMLBase`)*: Explicitly opts-out of DAE initialization.

### OrdinaryDiffEq.jl (Implicit ODE)

These methods from OrdinaryDiffEq are for `DAEProblem` specifications.

!!! note "v8: import from `OrdinaryDiffEqBDF`"

    `DImplicitEuler`, `DABDF2`, and `DFBDF` live in
    `OrdinaryDiffEqBDF` and are not in `OrdinaryDiffEq`'s default re-export
    set under v7. Load them with `using OrdinaryDiffEqBDF`.

  - `DImplicitEuler` *(from `OrdinaryDiffEqBDF`)* - 1st order A-L and stiffly stable adaptive implicit Euler
  - `DABDF2` *(from `OrdinaryDiffEqBDF`)* - 2nd order A-L stable adaptive BDF method.
  - `DFBDF` *(from `OrdinaryDiffEqBDF`)* - A fixed-leading coefficient adaptive-order adaptive-time BDF method,
    similar to `ode15i` or `IDA` in divided differences form.

### OrdinaryDiffEq.jl (Mass Matrix)

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

!!! note "v8: sublib mapping"

    Under OrdinaryDiffEq v7 only `Rosenbrock23`, `Rodas5P`, and `FBDF` from the
    list below are re-exported by the umbrella `OrdinaryDiffEq` package. The
    rest must be imported from their host sublib:

    | Section                    | Sublibrary                  |
    |----------------------------|-----------------------------|
    | Rosenbrock / Rosenbrock-W  | `OrdinaryDiffEqRosenbrock`  |
    | FIRK (`RadauIIA5`)         | `OrdinaryDiffEqFIRK`        |
    | SDIRK (`ImplicitEuler`, `ImplicitMidpoint`, `Trapezoid`) | `OrdinaryDiffEqSDIRK` |
    | Multistep (`QNDF`, `FBDF`, ...) | `OrdinaryDiffEqBDF`     |

!!! note
    
    The standard Hermite interpolation used for ODE methods in OrdinaryDiffEq.jl
    falls back to a linear interpolation on the differential variables. If the
    mass matrix is non-diagonal, the Hermite interpolation does not have a fallback
    and will error.

#### Rosenbrock Methods

  - `ROS3P` *(from `OrdinaryDiffEqRosenbrock`)* - 3rd order A-stable and stiffly stable Rosenbrock method. Keeps high
    accuracy on discretizations of nonlinear parabolic PDEs.
  - `Rodas3` *(from `OrdinaryDiffEqRosenbrock`)* - 3rd order A-stable and stiffly stable Rosenbrock method.
  - `RosShamp4` *(from `OrdinaryDiffEqRosenbrock`)*- An A-stable 4th order Rosenbrock method.
  - `Veldd4` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order D-stable Rosenbrock method.
  - `Velds4` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order A-stable Rosenbrock method.
  - `GRK4T` *(from `OrdinaryDiffEqRosenbrock`)* - An efficient 4th order Rosenbrock method.
  - `GRK4A` *(from `OrdinaryDiffEqRosenbrock`)* - An A-stable 4th order Rosenbrock method. Essentially "anti-L-stable"
    but efficient.
  - `Ros4LStab` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order L-stable Rosenbrock method.
  - `Rodas4` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order A-stable stiffly stable Rosenbrock method with a
    stiff-aware 3rd order interpolant
  - `Rodas42` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order A-stable stiffly stable Rosenbrock method with a
    stiff-aware 3rd order interpolant
  - `Rodas4P` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order A-stable stiffly stable Rosenbrock method with a
    stiff-aware 3rd order interpolant. 4th order on linear parabolic problems
    and 3rd order accurate on nonlinear parabolic problems (as opposed to lower
    if not corrected).
  - `Rodas4P2` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order L-stable stiffly stable Rosenbrock method with a stiff-aware
    3rd order interpolant. 4th order on linear parabolic problems and 3rd order accurate
    on nonlinear parabolic problems. It is an improvement of Roadas4P and in case of
    inexact Jacobians a second order W method.
  - `Rodas5` *(from `OrdinaryDiffEqRosenbrock`)* - A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware
    4th order interpolant.
  - `Rodas5P` *(from `OrdinaryDiffEqRosenbrock`)* - A 5th order A-stable stiffly stable Rosenbrock method with a stiff-aware
    4th order interpolant. Has improved stability in the adaptive time stepping embedding.

#### Rosenbrock-W Methods

  - `Rosenbrock23` *(from `OrdinaryDiffEqRosenbrock`)* - An Order 2/3 L-Stable Rosenbrock-W method which is good for very stiff equations with oscillations at low tolerances. 2nd order stiff-aware interpolation.
  - `Rosenbrock32` *(from `OrdinaryDiffEqRosenbrock`)* - An Order 3/2 A-Stable Rosenbrock-W method which is good for mildly stiff equations without oscillations at low tolerances. Note that this method is prone to instability in the presence of oscillations, so use with caution. 2nd order stiff-aware interpolation.
  - `RosenbrockW6S4OS` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order L-stable Rosenbrock-W method (fixed step only).
  - `ROS34PW1a` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order L-stable Rosenbrock-W method.
  - `ROS34PW1b` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order L-stable Rosenbrock-W method.
  - `ROS34PW2` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order stiffy accurate Rosenbrock-W method for PDAEs.
  - `ROS34PW3` *(from `OrdinaryDiffEqRosenbrock`)* - A 4th order strongly A-stable (Rinf~0.63) Rosenbrock-W method.

!!! note
    
    `Rosenbrock23` and `Rosenbrock32` have a stiff-aware interpolation but this interpolation is not safe for the algebraic variables.
    Thus use the interpolation (and therefore `saveat`) with caution if the default Hermite interpolation is used.

#### FIRK Methods

  - `RadauIIA5` *(from `OrdinaryDiffEqFIRK`)* - An A-B-L stable fully implicit Runge-Kutta method with internal
    tableau complex basis transform for efficiency.

#### SDIRK Methods

  - `ImplicitEuler` *(from `OrdinaryDiffEqSDIRK`)* - Stage order 1. A-B-L-stable. Adaptive
    timestepping through a divided differences estimate via memory. Strong-stability
    preserving (SSP).
  - `ImplicitMidpoint` *(from `OrdinaryDiffEqSDIRK`)* - Stage order 1. Symplectic. Good for when symplectic
    integration is required.
  - `Trapezoid` *(from `OrdinaryDiffEqSDIRK`)* - A second order A-stable symmetric ESDIRK method. "Almost
    symplectic" without numerical dampening. Also known as Crank-Nicolson when
    applied to PDEs. Adaptive timestepping via divided differences on the memory.
    Good for highly stiff equations which are non-oscillatory.

#### Multistep Methods

Quasi-constant stepping is the time stepping strategy which matches the classic
GEAR, LSODE,  and `ode15s` integrators. The variable-coefficient methods match
the ideas of the classic EPISODE integrator and early VODE designs. The Fixed
Leading Coefficient (FLC) methods match the behavior of the classic VODE and
Sundials CVODE integrator.

  - `QNDF1` *(from `OrdinaryDiffEqBDF`)* - An adaptive order 1 quasi-constant timestep L-stable numerical
    differentiation function (NDF) method. Optional parameter `kappa` defaults
    to Shampine's accuracy-optimal `-0.1850`.
  - `QBDF1` *(from `OrdinaryDiffEqBDF`)* - An adaptive order 1 L-stable BDF method. This is equivalent to
    implicit Euler but using the BDF error estimator.
  - `ABDF2` *(from `OrdinaryDiffEqBDF`)* - An adaptive order 2 L-stable fixed leading coefficient multistep
    BDF method.
  - `QNDF2` *(from `OrdinaryDiffEqBDF`)* - An adaptive order 2 quasi-constant timestep L-stable numerical
    differentiation function (NDF) method.
  - `QBDF2` *(from `OrdinaryDiffEqBDF`)* - An adaptive order 2 L-stable BDF method using quasi-constant timesteps.
  - `QNDF` *(from `OrdinaryDiffEqBDF`)* - An adaptive order quasi-constant timestep NDF method. Utilizes
    Shampine's accuracy-optimal `kappa` values as defaults (has a keyword argument
    for a tuple of `kappa` coefficients).
  - `QBDF` *(from `OrdinaryDiffEqBDF`)* - An adaptive order quasi-constant timestep BDF method.
  - `FBDF` *(from `OrdinaryDiffEqBDF`)* - A fixed-leading coefficient adaptive-order adaptive-time BDF method,
    similar to `ode15i` or `CVODE_BDF` in divided differences form.

### [Sundials.jl](@id dae_solve_sundials)

Note that this setup is not automatically included with DifferentialEquations.jl.
To use the following algorithms, you must install and use Sundials.jl:

```julia
using Pkg
Pkg.add("Sundials")
import Sundials
```

  - `IDA` *(from `Sundials`)*: A fixed-leading coefficient fully implicit BDF method. Efficient for large systems.

For more details on controlling the Sundials.jl solvers, see the
[Sundials detailed solver API page](@ref sundials)

### DASKR.jl

DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
using Pkg
Pkg.add("DASKR")
import DASKR
```

  - `daskr` *(from `DASKR`)* - This is a wrapper for the well-known DASKR algorithm.

For more details on controlling the DASKR.jl solvers, see the
[DASKR detailed solver API page](@ref daskr)

### DASSL.jl

  - `dassl` *(from `DASSL`)* - A native Julia implementation of the DASSL algorithm.

### ODEInterfaceDiffEq.jl

These methods require the DAE to be an `ODEProblem` in mass matrix form. For
extra options for the solvers, see the ODE solver page.

  - `seulex` *(from `ODEInterfaceDiffEq`)* - Extrapolation-algorithm based on the linear implicit Euler method.
  - `radau` *(from `ODEInterfaceDiffEq`)* - Implicit Runge-Kutta (Radau IIA) of variable order between 5 and 13.
  - `radau5` *(from `ODEInterfaceDiffEq`)* - Implicit Runge-Kutta method (Radau IIA) of order 5.
  - `rodas` *(from `ODEInterfaceDiffEq`)* - Rosenbrock 4(3) method.
