# [Financial Models](@id financial_models)

The financial models functionality is provided by DiffEqFinancial.jl and helps
the user build and solve the differential equation based financial models.

## SDE Model Library

The following constructors create `SDEProblem` types which can be solved using
the stochastic differential equation solvers.

### HestonProblem

```math
dS = μSdt + \sqrt{v}SdW_1 \\
dv = κ(Θ-v)dt + σ\sqrt{v}dW_2 \\
dW_1 dW_2 = ρ dt
```

Constructor:

```julia
HestonProblem(μ,κ,Θ,σ,ρ,u0,tspan)
```

### GeneralizedBlackScholesProblem

```math
d \ln S(t) = \left(r(t) - q(t) - \frac{Θ(t,S)^2}{2}\right)dt + σ dW_t
```

Solves for ``log S(t)``. Constructor:

```julia
GeneralizedBlackScholesProblem(r,q,Θ,σ,u0,tspan)
```

### BlackScholesProblem

```math
d \ln S(t) = \left(r(t) - \frac{Θ(t,S)^2}{2}\left)dt + σ dW_t
```

Solves for ``log S(t)``. Constructor:

```julia
BlackScholesProblem(r,Θ,σ,u0,tspan)
```

### ExtendedOrnsteinUhlenbeckProblem

```math
dx = a(b(t)-x)dt + σ dW_t
```

Constructor:

```julia
ExtendedOrnsteinUhlenbeckProblem(a,b,σ,u0,tspan)
```

### OrnsteinUhlenbeckProblem

```math
dx = a(r-x)dt + σ dW_t
```

Constructor:

```julia
OrnsteinUhlenbeckProblem(a,r,σ,u0,tspan)
```
