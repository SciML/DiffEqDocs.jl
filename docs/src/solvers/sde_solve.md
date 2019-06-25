# SDE Solvers

## Recommended Methods

For most Ito diagonal and scalar noise problems where a good amount of accuracy is
required and mild stiffness may be an issue, the `SOSRI` algorithm should
do well. If the problem has additive noise, then `SOSRA` will be the
optimal algorithm. At low tolerances (`<1e-4`?) `SRA3` will be more efficient,
though `SOSRA` is more robust to stiffness. For commutative noise, `RKMilCommute`
is a strong order 1.0 method which utilizes the commutivity property to greatly
speed up the Wiktorsson approximation and can choose between Ito and Stratonovich.
For non-commutative noise, difficult problems usually require adaptive time
stepping in order to be efficient. In this case, `LambaEM` and `LambaEulerHeun`
are adaptive and handle general non-diagonal problems (for Ito and Stratonovich
interpretations respectively). If adaptivity isn't necessary, the `EM` and
`EulerHeun` are good choices (for Ito and Stratonovich interpretations
respectively).

For stiff problems with additive noise, the high order adaptive method
`SKenCarp` is highly preferred and will solve problems with similar efficiency
as ODEs. If possible, stiff problems should be converted to make use of this
additive noise solver. If the noise term is large/stiff, then the split-step
methods are required in order for the implicit methods to be stable. For
Ito in this case, use `ISSEM` and for Stratonovich use `ISSEulerHeun`. These
two methods can handle any noise form.

If the noise term is not too large, for stiff problems with diagonal noise,
`ImplicitRKMil` is the most efficient method and can choose between Ito and
Stratonovich. For each of the theta methods, the parameter `theta` can be chosen.
The default is `theta=1/2` which will not dampen numerical oscillations
and thus is symmetric (and almost symplectic) and will lead to less error when
noise is sufficiently small. However, `theta=1/2` is not L-stable in the drift
term, and thus one can receive more stability (L-stability in the drift term)
with `theta=1`, but with a tradeoff of error efficiency in the low noise case.
In addition, the option `symplectic=true` will turns these methods into an
implicit Midpoint extension which is symplectic in distribution but has an
accuracy tradeoff.

## Mass Matrices and Stochastic DAEs

The stiff methods can solve stochastic equations with mass matrices
(including stochastic DAEs written in mass matrix form) when either
`symplectic=true` or `theta=1`. These methods interpret the mass matrix
equation as:

```math
Mu' = f(t,u)dt + Mg(t,u)dW_t
```

i.e. with no mass matrix inversion applied to the `g` term. Thus these methods
apply noise per dependent variable instead of on the combinations of the
dependent variables and this is designed for phenomenological noise on the
dependent variables (like multiplicative or additive noise)

## Special Noise Forms

Some solvers are for specialized forms of noise. Diagonal noise is the default
setup. Non-diagonal noise is specified via setting `noise_rate_prototype` to
a matrix in the `SDEProblem` type. A special form of non-diagonal noise,
commutative noise, occurs when the noise satisfies the following condition:

```math
\sum_{i=1}^d g_{i,j_1}(t,x) \frac{\partial g_{k,j_2}(t,x)}{\partial x_i} = \sum_{i=1}^d g_{i,j_2}(t,x) \frac{\partial g_{k,j_1}(t,x)}{\partial x_i}
```

for every ``j_1,j_2`` and ``k``. Additive noise is when ``g(t,u)=g(t)``,
i.e. is independent of `u`. Multiplicative noise is ``g_i(t,u)=a_i u``.

## Special Keyword Arguments

* `save_noise`: Determines whether the values of `W` are saved whenever the timeseries
  is saved. Defaults to true.
* `delta`: The `delta` adaptivity parameter for the natural error estimator.
  Determines the balance between drift and diffusion error. For more details, see
  [the publication](http://chrisrackauckas.com/assets/Papers/ChrisRackauckas-AdaptiveSRK.pdf).
* `seed`: Sets the seed for the random number generator. This overrides any seed
  set in the `SDEProblem`.

## Full List of Methods

### StochasticDiffEq.jl

Each of the StochasticDiffEq.jl solvers come with a linear interpolation.
Orders are given in terms of strong order.

#### Nonstiff Methods

- `EM`- The Euler-Maruyama method. Strong Order 0.5 in the Ito sense. Has an
  optional argument `split=true` for controlling step splitting. When splitting
  is enabled, the stability with large diffusion eigenvalues is improved. Can handle
  all forms of noise, including non-diagonal, scalar, and colored noise. Fixed
  time step only.†
- `LambaEM`- A modified Euler-Maruyama method with adaptive time stepping with
  an error estimator based on Lamba and Rackauckas. Has an optional argument
  `split=true` for controlling step splitting. When splitting is enabled, the
  stability with   large diffusion eigenvalues is improved. Strong Order 0.5 in
  the Ito sense. Can handle all forms of noise, including non-diagonal, scalar,
  and colored noise.†
- `EulerHeun` - The Euler-Heun method. Strong Order 0.5 in the Stratonovich sense.
  Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
  Fixed time step only.†
- `LambaEulerHeun` - A modified Euler-Heun method with adaptive time stepping
  with an error estimator based on Lamba due to Rackauckas. Strong order 0.5 in
  the Stratonovich sense. Can handle all forms of noise, including non-diagonal,
  scalar, and colored noise.†
- `RKMil` - An explicit Runge-Kutta discretization of the strong order 1.0
  Milstein method. Defaults to solving the Ito problem, but
  `RKMil(interpretation=:Stratonovich)` makes it solve the Stratonovich problem.
  Only handles scalar and diagonal noise.†
- `RKMilCommute` - An explicit Runge-Kutta discretization of the strong order 1.0
  Milstein method for commutative noise problems. Defaults to solving the Ito
  problem, but `RKMilCommute(interpretation=:Stratonovich)` makes it solve the
  Stratonovich problem. Uses a 1.5/2.0 error estimate for adaptive time stepping.†
- `WangLi3SMil_A` - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
- `WangLi3SMil_B` - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
- `WangLi3SMil_C` - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
- `WangLi3SMil_D` - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
- `WangLi3SMil_E` - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
- `WangLi3SMil_F` - fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
- `SRA` - Adaptive strong order 1.5 methods for additive Ito and Stratonovich SDEs.
  Default tableau is for SRA1. Can handle diagonal, non-diagonal and scalar
  additive noise.
- `SRI` - Adaptive strong order 1.5 methods for diagonal/scalar Ito SDEs.
  Default tableau is for SRIW1.
- `SRIW1` - Adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs.†
- `SRIW2` - Adaptive strong order 1.5 and weak order 3.0 for diagonal/scalar Ito SDEs.†
- `SOSRI` - Stability-optimized adaptive strong order 1.5 and weak order 2.0 for
  diagonal/scalar Ito SDEs. Stable at high tolerances and robust to stiffness.†
- `SOSRI2` - Stability-optimized adaptive strong order 1.5 and weak order 2.0 for
  diagonal/scalar Ito SDEs. Stable at high tolerances and robust to stiffness.†
- `SRA1` - Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak
  order 2. Can handle diagonal, non-diagonal, and scalar additive noise.†
- `SRA2` - Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak
  order 2. Can handle diagonal, non-diagonal, and scalar additive noise.†
- `SRA3` - Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak
  order 3. Can handle non-diagonal and scalar additive noise.†
- `SOSRA` - A stability-optimized adaptive SRA. Strong order 1.5 for additive Ito and
  Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal, and scalar
  additive noise. Stable at high tolerances and robust to stiffness.†
- `SOSRA2` - A stability-optimized adaptive SRA. Strong order 1.5 for additive Ito and
  Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal, and scalar
  additive noise. Stable at high tolerances and robust to stiffness.†

Example usage:

```julia
sol = solve(prob,SRIW1())
```
3-stage Milstein Methods `WangLi3SMil_A`, `WangLi3SMil_B`, `WangLi3SMil_D`, `WangLi3SMil_E` and `WangLi3SMil_F` are currently implemented for 1-Dimensional noise only.

#### Tableau Controls

For `SRA` and `SRI`, the following option is allowed:

* `tableau`: The tableau for an `:SRA` or `:SRI` algorithm. Defaults to SRIW1 or SRA1.

#### S-ROCK Methods

- `SROCK1` - is a fixed step size stabilized explicit method for stiff problems. Defaults to
  solving th Ito problem but `SROCK1(interpretation=:Stratonovich)` can make it solve
  the Stratonovich problem. Strong order of convergence is 0.5 and weak order 1, but is
  optimised to get order 1 in case os scalar/diagonal noise.
- `SROCKEM` - is fixed step Euler-Mayurama with first order ROCK stabilization thus can
  handle stiff problems. Only for Ito problems. Defaults to strong and weak order 1.0,
  but can solve with weak order 0.5 as `SROCKEM(strong_order_1=false)`.
- `SROCK2` - is a weak second order and strong first order fixed step stabilized method for
  stiff Ito problems.
- `SKSROCK` - is fixed step stabilized explicit method for stiff Ito problems. Strong order 0.5
  and weak order 1. This method has a better stability domain then `SROCK1`. Also it allows
  special post-processing techniques in case of ergodic dynamical systems, in the context of
  ergodic Brownian dynamics, to achieve order 2 accuracy.
- `TangXiaoSROCK2` - is a fixed step size stabilized expicit method for stiff problems. Only for
  Ito problems. Weak order of 2 and strog order of 1. Has 5 versions with different stability
  domains which can be used as `TangXiaoSROCK2(version_num=i)` where `i` is 1-5.

Note: Post processing in `SKSROCK` is still under development. Also `TangXiaoSROCK2` is also not fully
functional.

#### Stiff Methods

- `ImplicitEM` - An order 0.5 Ito drift-implicit method. This is a theta method which
  defaults to `theta=1/2` or the Trapezoid method on the drift term. This method
  defaults to `symplectic=false`, but when true and `theta=1/2` this is the
  implicit Midpoint method on the drift term and is symplectic in distribution.
  Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
  Uses a 1.0/1.5 heuristic for adaptive time stepping.
- `ImplicitEulerHeun` - An order 0.5 Stratonovich drift-implicit method. This is a
  theta method which defaults to `theta=1/2` or the Trapezoid method on the
  drift term. This method defaults to `symplectic=false`, but when true and
  `theta=1/2` this is the implicit Midpoint method on the drift term and is
  symplectic in distribution. Can handle all forms of noise, including
  non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for
  adaptive time stepping.
- `ImplicitRKMil` - An order 1.0 drift-implicit method. This is a theta method which
  defaults to `theta=1/2` or the Trapezoid method on the drift term. Defaults
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
  distribution. Can handle all forms of noise, including non-diagonal,Q scalar,
  and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.
- `SKenCarp` - Adaptive L-stable drfit-implicit strong order 1.5 for additive
  Ito and Stratonovich SDEs with weak order 2. Can handle diagonal, non-diagonal
  and scalar additive noise.†

#### Derivative-Based Methods

The following methods require analytic derivatives of the diffusion term.

- `PCEuler` - The predictor corrector euler method. Strong Order 0.5 in the Ito
  sense. Requires the ggprime function, which is defined as
  ```math
    \text{ggprime}^k(t,x) = \sum_{j=1}^m \sum_{i=1}^d g_{i,j}(t,x) \frac{\partial g_{k,j}(t,x)}{\partial x_i}.
  ```
  This can also be understood more intuitively in vector/matrix form as,
  ```math
  \text{ggprime}(t,x) = \sum_{j=1}^m \bar{\mathcal{J}}[\vec g^{(j)}(t,x)] \vec g^{(j)}(t,x).
  ```
  where ``\vec g^{(j)}`` is the noise vector for the j'th noise channel and ``\bar{\mathcal{J}}`` is the Jacobian of the j'th   noise vector.

  The default settings for the drift implicitness is `theta=0.5` and
  the diffusion implicitness is `eta=0.5`.  

### StochasticCompositeAlgorithm

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

### SimpleDiffEq.jl

This setup provides access to simplified versions of a few SDE solvers. They
mostly exist for experimentation, but offer shorter compile times. They have
limitations compared to StochasticDiffEq.jl.

  - `SimpleEM` - A fixed timestep solve method for Euler-Maruyama. Only works
    with non-colored Gaussian noise.

Note that this setup is not automatically included with DifferentialEquaitons.jl.
To use the following algorithms, you must install and use SimpleDiffEq.jl:

```julia
]add SimpleDiffEq
using SimpleDiffEq
```

### BridgeDiffEq.jl

Bridge.jl is a set of fixed timestep algorithms written in Julia. These methods
are made and optimized for out-of-place functions on immutable (static vector)
types. Note that this setup is not automatically included with
DifferentialEquaitons.jl. To use the following algorithms, you must install and
use BridgeDiffEq.jl:

```julia
Pkg.clone("https://github.com/JuliaDiffEq/BridgeDiffEq.jl")
using BridgeDiffEq
```

- `BridgeEuler` - Strong order 0.5 Euler-Maruyama method for Ito equations.†
- `BridgeHeun` - Strong order 0.5 Euler-Heun method for Stratonovich equations.†
- `BridgeSRK` - Strong order 1.0 derivative-free stochastic Runge-Kutta method
  for scalar (`<:Number`) Ito equations.†

##### Notes

†: Does not step to the interval endpoint. This can cause issues with discontinuity
detection, and [discrete variables need to be updated appropriately](../../features/diffeq_arrays.html).
