# DAE Solvers

## Recomended Methods

Currently, the only method in the ecosystem is `IDA`.

## Special Keyword Arguments

## Implemented Solvers

### Sundials.jl

- `IDA` - This is the IDA method from the Sundials.jl package.

### DASKR.jl

DASKR.jl is not automatically included by DifferentialEquations.jl. To use this
algorithm, you will need to install and use the package:

```julia
Pkg.add("DASKR")
using DASKR
```

- `DASRK` - This is a wrapper for the well-known DASKR algorithm.
