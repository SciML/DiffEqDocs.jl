# [DE types](@id de_types)


## Discrete

```@docs
DiffEqBase.AbstractDiscreteFunction
DiscreteFunction
DiffEqBase.AbstractDiscreteProblem
DiscreteProblem
```

## ODE

```@docs
DiffEqBase.AbstractODEFunction
ODEFunction
DiffEqBase.AbstractODEProblem
ODEProblem
DiffEqBase.StandardODEProblem
DiffEqBase.AbstractODESolution
DiffEqBase.ODESolution
DiffEqBase.AbstractODEAlgorithm
DiffEqBase.AbstractODEIntegrator
```

### Dynamical ODEs

```@docs
DynamicalODEFunction
DiffEqBase.AbstractDynamicalODEProblem
DynamicalODEProblem
```

### Second-order ODEs

```@docs
DiffEqBase.AbstractSecondOrderODEProblem
SecondOrderODEProblem
DiffEqBase.AbstractSecondOrderODEAlgorithm
DiffEqBase.AbstractSecondOrderODEIntegrator
```

### Split ODEs

```@docs
SplitFunction
DiffEqBase.AbstractSplitODEProblem
SplitODEProblem
```

### Steady state problems

```@docs
DiffEqBase.AbstractSteadyStateProblem
SteadyStateProblem
DiffEqBase.AbstractSteadyStateSolution
DiffEqBase.SteadyStateSolution
DiffEqBase.AbstractSteadyStateAlgorithm
DiffEqBase.AbstractSteadyStateIntegrator
```

### Boundary value problems

```@docs
DiffEqBase.TwoPointBVPFunction
DiffEqBase.AbstractBVProblem
DiffEqBase.StandardBVProblem
BVProblem
TwoPointBVProblem
```

### Analytical problems

```@docs
AbstractAnalyticalProblem
AnalyticalProblem
DiffEqBase.AbstractAnalyticalSolution
```


## SDE

SDE problems are subtypes of RODE problems.

```@docs
DiffEqBase.AbstractSDEFunction
SDEFunction
DiffEqBase.AbstractSDEProblem
DiffEqBase.StandardSDEProblem
SDEProblem
DiffEqBase.AbstractSDEAlgorithm
DiffEqBase.AbstractSDEIntegrator
```

### Split SDEs

```@docs
SplitSDEFunction
DiffEqBase.AbstractSplitSDEProblem
SplitSDEProblem
```


## RODE

```@docs
DiffEqBase.AbstractRODEFunction
RODEFunction
DiffEqBase.AbstractRODEProblem
RODEProblem
DiffEqBase.AbstractRODESolution
DiffEqBase.RODESolution
DiffEqBase.AbstractRODEAlgorithm
DiffEqBase.AbstractRODEIntegrator
```


## DDE

```@docs
DiffEqBase.AbstractDDEFunction
DDEFunction
DiffEqBase.AbstractDDEProblem
DDEProblem
DiffEqBase.AbstractConstantLagDDEProblem
DiffEqBase.AbstractDDESolution
DiffEqBase.AbstractDDEAlgorithm
DiffEqBase.AbstractDDEIntegrator
DiffEqBase.AbstractHistoryFunction
```


## SDDE

```@docs
DiffEqBase.AbstractSDDEFunction
SDDEFunction
DiffEqBase.AbstractSDDEProblem
SDDEProblem
DiffEqBase.AbstractConstantLagSDDEProblem
DiffEqBase.AbstractSDDEAlgorithm
DiffEqBase.AbstractSDDEIntegrator
```


## DAE

```@docs
DiffEqBase.AbstractDAEFunction
DAEFunction
DiffEqBase.AbstractDAEProblem
DAEProblem
DiffEqBase.AbstractDAESolution
DiffEqBase.DAESolution
DiffEqBase.AbstractDAEAlgorithm
DiffEqBase.AbstractDAEIntegrator
```


### PDE

```@docs
DiffEqBase.AbstractPDEProblem
PDEProblem
```


## Jump problems

```@docs
DiffEqBase.AbstractJumpProblem
```


## Noise problems

```@docs
DiffEqBase.AbstractNoiseProblem
NoiseProblem
```


## Basic problems (for testing?)

### Linear

```@docs
DiffEqBase.AbstractLinearProblem
DiffEqBase.LinearProblem
DiffEqBase.AbstractLinearSolution
DiffEqBase.LinearSolution
DiffEqBase.AbstractLinearAlgorithm
```

### Nonlinear

```@docs
DiffEqBase.AbstractNonlinearProblem
DiffEqBase.NonlinearProblem
DiffEqBase.AbstractNonlinearSolution
DiffEqBase.AbstractNonlinearAlgorithm
```

### Quadrature

```@docs
DiffEqBase.AbstractQuadratureProblem
DiffEqBase.QuadratureProblem
DiffEqBase.AbstractQuadratureSolution
DiffEqBase.AbstractQuadratureAlgorithm
```


## Sensitivity problems
```@docs
DiffEqBase.DESensitivity
DiffEqBase.AbstractSensitivitySolution
```
