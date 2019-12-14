# Callbacks

```@meta
CurrentModule = DiffEqBase
```

## Types

```@docs
DiffEqBase.DECallback
DiffEqBase.AbstractDiscreteCallback
DiffEqBase.AbstractContinuousCallback
DiffEqBase.CallbackCache
```


## Functions

```@docs
DiffEqBase.split_callbacks
DiffEqBase.initialize!(cb::CallbackSet, u, t, integrator::DEIntegrator)
```
