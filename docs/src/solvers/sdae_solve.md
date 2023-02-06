# [SDAE Solvers](@id sdae_solve)

## Recommended Methods

The recommendations for SDAEs are the same recommended implicit SDE methods for
stiff equations when the SDAE is specified in mass matrix form.

#### Mass Matrix Form

  - `ImplicitEM` - An order 0.5 Ito drift-implicit method. This is a theta method which
    defaults to `theta=1` or the Trapezoid method on the drift term. This method
    defaults to `symplectic=false`, but when true and `theta=1/2` this is the
    implicit Midpoint method on the drift term and is symplectic in distribution.
    Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
    Uses a 1.0/1.5 heuristic for adaptive time stepping.
  - `STrapezoid` - An alias for `ImplicitEM` with `theta=1/2`
  - `SImplicitMidpoint` - An alias for `ImplicitEM` with `theta=1/2` and `symplectic=true`
  - `ImplicitEulerHeun` - An order 0.5 Stratonovich drift-implicit method. This is a
    theta method which defaults to `theta=1/2` or the Trapezoid method on the
    drift term. This method defaults to `symplectic=false`, but when true and
    `theta=1` this is the implicit Midpoint method on the drift term and is
    symplectic in distribution. Can handle all forms of noise, including
    non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for
    adaptive time stepping.
  - `ImplicitRKMil` - An order 1.0 drift-implicit method. This is a theta method which
    defaults to `theta=1` or the Trapezoid method on the drift term. Defaults
    to solving the Ito problem, but `ImplicitRKMil(interpretation=:Stratonovich)`
    makes it solve the Stratonovich problem. This method defaults to
    `symplectic=false`, but when true and `theta=1/2` this is the
    implicit Midpoint method on the drift term and is symplectic in distribution.
    Handles diagonal and scalar noise. Uses a 1.5/2.0 heuristic for adaptive
    time stepping.
  - `ISSEM` - An order 0.5 split-step Ito implicit method. It is fully implicit,
    meaning it can handle stiffness in the noise term. This is a theta method which
    defaults to `theta=1` or the Trapezoid method on the drift term. This method
    defaults to `symplectic=false`, but when true and `theta=1/2` this is the
    implicit Midpoint method on the drift term and is symplectic in distribution.
    Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
    Uses a 1.0/1.5 heuristic for adaptive time stepping.
  - `ISSEulerHeun` - An order 0.5 split-step Stratonovich implicit method. It is
    fully implicit, meaning it can handle stiffness in the noise term. This is a
    theta method which defaults to `theta=1` or the Trapezoid method on the drift
    term. This method defaults to `symplectic=false`, but when true and `theta=1/2`
    this is the implicit Midpoint method on the drift term and is symplectic in
    distribution. Can handle all forms of noise, including non-diagonal, Q scalar,
    and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.
  - `SKenCarp` - Adaptive L-stable drift-implicit strong order 1.5 for additive
    Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal
    and scalar additive noise.\*†

## Notes

†: Does not step to the interval endpoint. This can cause issues with discontinuity
detection, and [discrete variables need to be updated appropriately](@ref diffeq_arrays).

\*:  Note that although `SKenCarp` uses the same table as `KenCarp3`, solving a ODE problem using `SKenCarp` by setting `g(du,u,p,t) = du .= 0` will take many more steps than `KenCarp3` because error estimator of `SKenCarp` is different (because of noise terms) and default value of `qmax` (maximum permissible ratio of relaxing/tightening `dt` for adaptive steps) is smaller for StochasticDiffEq algorithms.
