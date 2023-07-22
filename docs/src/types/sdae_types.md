# SDAE Problems

## Mathematical Specification of a Stochastic Differential-Algebraic Equation (SDAE) Problem

To define an SDAE, you simply define an SDE Problem with the forcing function `f`,
the noise function `g`, a mass matrix `M` and the initial condition `u₀` which
define the SDAE in mass matrix form:

```math
M du = f(u,p,t)dt + Σgᵢ(u,p,t)dWⁱ
```

`f` and `g` should be specified as `f(u,p,t)` and  `g(u,p,t)` respectively, and `u₀`
should be an AbstractArray whose geometry matches the desired geometry of `u`.
Note that we are not limited to numbers or vectors for `u₀`; one is allowed to
provide `u₀` as arbitrary matrices / higher dimension tensors as well. A vector
of `g`s can also be defined to determine an SDE of higher Ito dimension.

Nonsingular mass matrices correspond to constraint equations and thus a stochastic
DAE.

## Example

```julia
const mm_A = [-2.0 1 4
    4 -2 1
    0 0 0]
function f!(du, u, p, t)
    du[1] = u[1]
    du[2] = u[2]
    du[3] = u[1] + u[2] + u[3] - 1
end

function g!(du, u, p, t)
    @. du = 0.1
end

prob = SDEProblem(SDEFunction(f!, g!; mass_matrix = mm_A), g!,
    ones(3), (0.0, 1.0))
```
