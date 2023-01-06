# SDDE Solvers

`solve(prob::AbstractSDDEProblem, alg; kwargs)`

Solves the SDDE defined by `prob` using the algorithm `alg`. If no algorithm is
given, a default algorithm will be chosen.

## Recommended Methods

The recommended method for SDDE problems are the `SDE` algorithms. On SDEs you
simply reuse the same algorithm as the `SDE` solver, and StochasticDelayDiffEq.jl
will convert it to an SDDE solver. The recommendations for SDDE solvers match
those of SDEs, except that only up to strong order 1 is recommended. Also note
that order 1 is currently only attainable if there is no delay term in the
diffusion function ``g``: delays in the drift function ``f`` are compatible
with first order convergence. Theoretical issues with higher order methods
(1.5+) on SDDEs is currently unknown.

Note that adaptive time stepping utilizes the same
[rejection sampling with memory](https://chrisrackauckas.com/assets/Papers/ChrisRackauckas-AdaptiveSRK.pdf)
technique as SDEs, but no proof of convergence is known for SDDEs.

## Example

```julia
function hayes_modelf(du,u,h,p,t)
    τ,a,b,c,α,β,γ = p
    du .= a.*u .+ b .* h(p,t-τ) .+ c
end
function hayes_modelg(du,u,h,p,t)
    τ,a,b,c,α,β,γ = p
    du .= α.*u .+ γ
end
h(p,t) = (ones(1) .+ t);
tspan = (0.,10.)

pmul = [1.0,-4.,-2.,10.,-1.3,-1.2, 1.1]
padd = [1.0,-4.,-2.,10.,-0.0,-0.0, 0.1]

prob = SDDEProblem(hayes_modelf, hayes_modelg, [1.], h, tspan, pmul; constant_lags = (pmul[1],));
sol = solve(prob,RKMil())
```
