# SDE Solvers

## Implemented Solvers

In addition to the standard Euler-Maruyama method, specialized versions of higher
order Runge-Kutta methods are implemented which give increased accuracy and speed.

  * Euler-Maruyama
  * Milstein
  * Rossler-SRK

## Solver Documentation

`solve(prob::SDEProblem,tspan)`

Solves the SDE as defined by prob on the time interval tspan. If not given, tspan defaults to [0,1].

### Special Keyword Arguments

* `discard_length` - Size at which to discard future information in adaptive. Default is 1e-15.
* `tableau`: The tableau for an `:SRA` or `:SRI` algorithm. Defaults to SRIW1 or SRA1.
* `adaptivealg`: The adaptive timestepping algorithm. Default is `:RSwm3`.
* `alg`: String which defines the solver algorithm. Defult is "SRIW1Optimized". Possibilities are:

    - `:EM`- The Euler-Maruyama method.
    - `:RKMil` - An explicit Runge-Kutta discretization of the strong Order 1.0 Milstein method.
    - `:SRA` - The strong Order 2.0 methods for additive SDEs due to Rossler. Not yet implemented.
      Default tableau is for SRA1.
    - `:SRI` - The strong Order 1.5 methods for diagonal/scalar SDEs due to Rossler.
      Default tableau is for SRIW1.
    - `:SRIW1Optimized` - An optimized version of SRIW1. Strong Order 1.5.
    - `:SRA1Optimized` - An optimized version of SRIA1. Strong Order 2.0.
    - `:SRAVectorized` - A vectorized implementation of SRA algorithms. Requires 1-dimensional problem.
    - `:SRIVectorized` - A vectorized implementation of SRI algorithms. Requires 1-dimensional problem.
