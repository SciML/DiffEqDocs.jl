# SDE Solvers

## Recommended Methods

For most diagonal and scalar noise problems where a good amount of accuracy is
required and stiffness may be an issue, the `SRIW1Optimized` algorithm should do
well. If the problem has additive noise, then `SRA1Optimized` will be the optimal
algorithm. For non-commutative noise, `EM` and `EulerHeun` will be the most accurate
(for Ito and Stratonovich interpretations respectively).

## Special Keyword Arguments

* `save_noise`: Determines whether the values of `W` are saved whenever the timeseries
  is saved. Defaults to true.
* `delta`: The `delta` adaptivity parameter for the natural error estimator. For
  more details, see [the publication](http://chrisrackauckas.com/assets/Papers/ChrisRackauckas-AdaptiveSRK.pdf).

# Full List of Methods

## StochasticDiffEq.jl

Each of the StochasticDiffEq.jl solvers come with a linear interpolation.

- `EM`- The Euler-Maruyama method. Strong Order 0.5 in the Ito sense.†
- `EulerHeun` - The Euler-Heun method. Strong Order 0.5 in the Stratonovich sense.
- `RKMil` - An explicit Runge-Kutta discretization of the strong Order 1.0 (Ito) Milstein method.†
- `SRA` - The strong Order 2.0 methods for additive Ito and Stratonovich SDEs due to Rossler.
  Default tableau is for SRA1.
- `SRI` - The strong Order 1.5 methods for diagonal/scalar Ito SDEs due to Rossler.
  Default tableau is for SRIW1.
- `SRIW1` - An optimized version of SRIW1. Strong Order 1.5 for diagonal/scalar Ito SDEs.†
- `SRA1` - An optimized version of SRA1. Strong Order 2.0 for additive Ito and Stratonovich SDEs.†

Example usage:

```julia
sol = solve(prob,SRIW1())
```

For `SRA` and `SRI`, the following option is allowed:

* `tableau`: The tableau for an `:SRA` or `:SRI` algorithm. Defaults to SRIW1 or SRA1.

†: Does not step to the interval endpoint. This can cause issues with discontinuity
detection, and [discrete variables need to be updated appropriately](../features/diffeq_arrays.html).

#### StochasticCompositeAlgorithm

One unique feature of StochasticDiffEq.jl is the `StochasticCompositeAlgorithm`, which allows
you to, with very minimal overhead, design a multimethod which switches between
chosen algorithms as needed. The syntax is `StochasticCompositeAlgorithm(algtup,choice_function)`
where `algtup` is a tuple of StochasticDiffEq.jl algorithms, and `choice_function`
is a function which declares which method to use in the following step. For example,
we can design a multimethod which uses `EM()` but switches to `RKMil()` whenever
`dt` is too small:

```julia
choice_function(integrator) = (Int(integrator.dt<0.001) + 1)
alg_switch = StochasticCompositeAlgorithm((EM(),RKMil()),choice_function)
```

The `choice_function` takes in an `integrator` and thus all of the features
available in the [Integrator Interface](@ref)
can be used in the choice function.
