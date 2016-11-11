# Differential Algebraic Equation Solvers

`solve(prob::DAEProblem,tspan)`

Solves the DAE as defined by prob on the time interval tspan. If not given, tspan defaults to [0,1].

### Special Keyword Arguments

* `alg`: String which defines the solver algorithm. Default is "idasol". Possibilities are:
  - `idasol`: The DAE solver from Sundials
