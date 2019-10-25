# DiffEqBase API Overview


## DE types

`DiffEqBase.jl` defines specialized types and functions for several categories
of differential equation problems (ODEs, SDEs, etc.). Most categories have
definitions for the following:

* Function types
  (`<: `[`AbstractDiffEqFunction`](@ref DiffEqBase.AbstractDiffEqFunction)),
  which wrap a Julia function along with additional data such as a mass matrix
  or additional functions for calculating Jacobians.
* Problem types (`<:`[`DEProblem`](@ref DiffEqBase.DEProblem)), which represent
  a problem to be solved. They include a function of the corresponding type
  along with information such as initial conditions, parameter values, and the
  time span to solve.
* Algorithm types (`<: `[`DEAlgorithm`](@ref DiffEqBase.DEAlgorithm)), which solve
  `DEProblem`s.
* Solution types (`<: `[`DESolution`](@ref DiffEqBase.DESolution)), which are
  the result of an algorithm solving a `DEProblem`.

See [DE Types](@ref de_types) for documentation of specialized code for each problem type.


## Type parameters

The following are common type parameters that appear throughout the package:

* `uType` - Element type of the state vector of a DE function or problem.
  Typically a subtype of `Real`.
* `tType` - Type of time variable used in a DE function or problem. Typically a
  subtype of `Real`.
* `isinplace` - Boolean value which indicates whether a DE function operates
  in-place (sets values of an array passed as its first argument rather than
  returning a new array). Also appears as `iip`.
* `ND`


## Index

```@index
```
