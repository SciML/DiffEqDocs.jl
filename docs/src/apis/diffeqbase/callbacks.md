# Callbacks

```@meta
CurrentModule = DiffEqBase
```

## Types

```@docs
DiffEqBase.DECallback
DiffEqBase.AbstractDiscreteCallback
DiscreteCallback
DiffEqBase.AbstractContinuousCallback
ContinuousCallback
VectorContinuousCallback
CallbackSet
DiffEqBase.CallbackCache
```


## Functions

```@docs
DiffEqBase.split_callbacks
DiffEqBase.initialize!(cb::CallbackSet, u, t, integrator::DEIntegrator)
```
