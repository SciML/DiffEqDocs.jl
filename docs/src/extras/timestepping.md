# [Timestepping Method Descriptions](@id timestepping)

```@meta
CurrentModule = OrdinaryDiffEqCore
```

## Common Setup

All methods start by calculating a scaled error estimate on each scalar component of ``u``:

```math
err^{scaled}_i = norm(err_i/(abstol_i + max(uprev_i,u_i)reltol_i))
```

On this scaled error estimate, we calculate the norm. This norm is usually the
Hairer norm:

```math
norm(x) = sqrt(sum(x^2)/length(x))
```

This norm works well because it does not change if we add new pieces to the
differential equation: it scales our error by the number of equations so that
independent equations will not step differently than a single solve.

In all cases, the step is rejected if ``err^{scaled}>1`` since that means the
error is larger than the tolerances, and the step is accepted if
``err^{scaled}<1``.

## Integral Controller (Standard Controller)

The integral control algorithm is the “standard algorithm” for adaptive
timestepping. Note that it is not the default in DifferentialEquations.jl
because it is usually awful for performance, but it is explained first because
it is the most widely taught algorithm and others build on its techniques.

The control simply changes `dt` proportional to the error. There is an
exponentiation based on the order of the algorithm which goes back to a result
by Cechino for the optimal stepsize to reduce the error. The algorithm is:

```julia
qtmp = integrator.EEst^(1 / (alg_adaptive_order(integrator.alg) + 1)) /
       integrator.opts.gamma
@fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
integrator.dtnew = integrator.dt / q
```

Thus, `q` is the scaling factor for `dt`, and it must be between `qmin` and `qmax`.
`gamma` is the safety factor, `0.9`, for how much `dt` is decreased below the
theoretical “optimal” value.

Since proportional control is “jagged”, i.e. can cause large changes between
one step to the next, it can effect the stability of explicit methods. Thus,
it's only applied by default to low order implicit solvers.

```@docs
OrdinaryDiffEqCore.IController
```

## Proportional-Integral Controller (PI Controller)

The proportional-integral control algorithm is a standard control algorithm
from control theory. It mixes proportional control with memory in order to
make the timesteps more stable, which actually increases the adaptive stability
region of the algorithm. This stability property means that it's well-suited
for explicit solvers, and it's applied by default to the Rosenbrock methods
as well. The form for the updates is:

```julia
EEst, beta1, q11, qold, beta2 = integrator.EEst, integrator.opts.beta1, integrator.q11,
integrator.qold, integrator.opts.beta2
@fastmath q11 = EEst^beta1
@fastmath q = q11 / (qold^beta2)
integrator.q11 = q11
@fastmath q = max(inv(integrator.opts.qmax),
    min(inv(integrator.opts.qmin), q / integrator.opts.gamma))
if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
end
q
```

`beta1` is the gain on the proportional part, and `beta2` is the gain for the
history portion. `qoldinit` is the initialized value for the gain history.

```@docs
OrdinaryDiffEqCore.PIController
```

## Proportional-Integral-Derivative Controller (PID Controller)

```@docs
OrdinaryDiffEqCore.PIDController
```

## Gustafsson Acceleration

The Gustafsson acceleration algorithm accelerates changes so that way algorithms
can more swiftly change to handle quick transients. This algorithm is thus
well-suited for stiff solvers where this can be expected, and is the default
for algorithms like the (E)SDIRK methods.

```julia
gamma = integrator.opts.gamma
niters = integrator.cache.newton_iters
fac = min(gamma,
    (1 + 2 * integrator.alg.max_newton_iter) * gamma /
    (niters + 2 * integrator.alg.max_newton_iter))
expo = 1 / (alg_order(integrator.alg) + 1)
qtmp = (integrator.EEst^expo) / fac
@fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
end
integrator.qold = q
q
```

In this case, `niters` is the number of Newton iterations which was required in
the most recent step of the algorithm. Note that these values are used differently
depending on acceptance and rejection. When the step is accepted, the
following logic is applied:

```julia
if integrator.success_iter > 0
    expo = 1 / (alg_adaptive_order(integrator.alg) + 1)
    qgus = (integrator.dtacc / integrator.dt) *
           (((integrator.EEst^2) / integrator.erracc)^expo)
    qgus = max(inv(integrator.opts.qmax),
        min(inv(integrator.opts.qmin), qgus / integrator.opts.gamma))
    qacc = max(q, qgus)
else
    qacc = q
end
integrator.dtacc = integrator.dt
integrator.erracc = max(1e-2, integrator.EEst)
integrator.dt / qacc
```

When it rejects, it is the same as the proportional control:

```julia
if integrator.success_iter == 0
    integrator.dt *= 0.1
else
    integrator.dt = integrator.dt / integrator.qold
end
```

```@docs
OrdinaryDiffEqCore.PredictiveController
```

## Abstract Controller

The `AbstractController` type allows one to implement custom adaptive
timestepping schemes easily. This can be useful if one wants to use
a priori error estimates, for instance.

To implement a custom controller, subtype `AbstractController` for the
specific DE system.

```julia
struct CustomController <: AbstractController
end
```

and overload

```julia
function stepsize_controller!(integrator, controller::CustomController, alg)
    ...
    nothing
end
function step_accept_controller!(integrator, controller::CustomController, alg)
    ...
    nothing
end
function step_reject_controller!(integrator, controller::CustomController, alg)
    ...
    nothing
end
```

For instance, the PI controller for SDEs can be reproduced by

```julia
struct CustomController <: StochasticDiffEq.AbstractController
end

function StochasticDiffEq.stepsize_controller!(integrator::StochasticDiffEq.SDEIntegrator,
        controller::CustomController, alg)
    integrator.q11 = DiffEqBase.value(FastPower.fastpower(
        integrator.EEst, controller.beta1))
    integrator.q = DiffEqBase.value(integrator.q11 /
                                    FastPower.fastpower(integrator.qold, controller.beta2))
    integrator.q = DiffEqBase.value(max(inv(integrator.opts.qmax),
        min(inv(integrator.opts.qmin),
            integrator.q / integrator.opts.gamma)))
    nothing
end

function StochasticDiffEq.step_accept_controller!(
        integrator::StochasticDiffEq.SDEIntegrator,
        controller::CustomController, alg)
    integrator.dtnew = DiffEqBase.value(integrator.dt / integrator.q) *
                       oneunit(integrator.dt)
    nothing
end

function step_reject_controller!(integrator::StochasticDiffEq.SDEIntegrator,
        controller::CustomController, alg)
    integrator.dtnew = integrator.dt / min(inv(integrator.opts.qmin),
        integrator.q11 / integrator.opts.gamma)
end
```

and used via

```julia
# Assuming: import DifferentialEquations as DE
sol = DE.solve(prob, DE.EM(), dt = dt, adaptive = true, controller = CustomController())
```
