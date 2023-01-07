# Solver Compatibility Chart

This chart is for documenting the compatibility of the component solver packages
to the common interface. An `x` means that the option is implemented or the
add-on functionality will work with the given solver. A blank means that
the option has not been implemented, or that a given add-on has not been tested
with a given package. If there are any errors in this chart, please file an
issue or submit a pull-request.

| Option                                 | OrdinaryDiffEq.jl | Sundials.jl | ODE.jl | ODEInterface.jl | LSODA.jl | StochasticDiffEq.jl | DelayDiffEq.jl | DASKR.jl | DASSL.jl
|----------------------------------------|-------------------|-------------|--------|-----------------|----------|---------------------|----------------|----------|----------
| Nonlinear Dense (continuous) output    | x                 | x           |        |                 |          | x                   | x              | x        |          
| Tolerance control                      | x                 | x           | x      | x               | x        | x                   | x              | x        | x        
| Advanced stepsize control              | x                 | 0           |        | x               | 0        | x                   | x              | 0        |          
| Mass Matrices^                         | x                 | 0           |        | x               | 0        | x                   | x              | 0        |     
| Analytical Jacobians^†                 | x                 | x           |        | x               |          | x                   | x              | x        |     
| General Performance Overloads^†        | x                 | 0           |        | 0               | 0        | x                   | x              | 0        |  
| internalnorm                           | x                 | 0           | x      | 0               | 0        | x                   | x              | 0        |          
| Initial dt                             | x                 | x           | x      | x               |          | x                   | x              | x        |          
| save_everystep                         | x                 | x           | x      | x               | x        | x                   | x              | x        |          
| saveat                                 | x                 | x           | x      | x               | x        | x                   | x              | x        |          
| tstops                                 | x                 | x           |        | 0               |          | x                   | x              | x        |          
| d_discontinuities                      | x                 |             |        | 0               |          | x                   | x              |          |          
| isoutofdomain                          | x                 |             | x      |                 |          | x                   | x              |          |          
| Allows reverse time direction          | x                 | x           | x      | x               | x        | x                   | x              |          |          
| Unitful numbers                        | x                 | 0           |        | 0               | 0        |                     | x              | 0        |          
| Arbitrary dimension arrays             | x                 | x           | x      | x               | x        | x                   | x              | x        | x        
| Complex numbers                        | p                 |             |        |                 |          | x                   | p              |          |          
| Arbitrary precision                    | x                 | 0           | x      | 0               | 0        | x                   | x              | 0        | x        
| ApproxFun types                        | x                 | 0           |        | 0               | 0        |                     | x              | 0        |          
| Progress monitoring                    | x                 |             |        |                 |          | x                   | x              |          |          
| Integrator interface                   | x                 | x           |        | 0               |          | x                   | x              |          |          
| Resizability                           | x                 | 0           |        | 0               | 0        | x                   | x              | 0        |          
| Cache iterator                         | x                 | 0           |        | 0               | 0        | x                   | x              | 0        |          
| Can choose linear solvers              | x                 | s           |        |                 |          | x                   | x              | s        | x        
| Can choose nonlinear solvers           | x                 | 0           |        | 0               | 0        | x                   | x              | 0        | x        
| Can use out of place natively          | x                 | 0           | x      | 0               | 0        | x                   | x              | 0        | x        
| Can use inplace natively               | x                 | x           |        | x               | x        | x                   | x              | x        |         
| Compatible with DiffEqDevTools         | x                 | x           | x      | x               | x        | x                   | x              | x        |          
| Compatible with ParameterizedFunctions | x                 | x           | x      | x               | x        | x                   | x              | x        |          
| Continuous Callbacks                   | x                 | x           |        | x               |          | x                   | x              |          | x        
| Discrete Callbacks                     | x                 | x           |        | x               |          | x                   | x              |          |          
| Monte Carlo Simulations                | x                 | x           | x      | x               | x        | x                   | x              | x        |          
| Parameter Estimation                   | x                 | n           | n      | n               | n        | x                   | x              | n        | x        
| Parameter Sensitivity Analysis         | x                 | x           | x      | x               | x        |                     | x              |          |          
| Plotting and solution handling         | x                 | x           | x      | x               | x        | x                   | x              | x        | x          

* x: Full compatibility
* p: Partial compatibility, only in nonstiff methods, unless the Jacobian is provided.
* n: General compatibility, but not compatible with routines which
  require being able to autodifferentiate through the entire solver.
* 0: Not possible. This is generally due to underlying inflexibility in a wrapped
  library.
* s: Special, Sundials has its own linear solver choices.
* ^: Only stiff (implicit) methods.
* †: For packages with compatibility, no warning is given when a specific algorithm
  does not need to use this feature.

All blank spaces are possible future additions.
