# Ensembles


## Types

```@docs
DiffEqBase.AbstractEnsembleProblem
EnsembleProblem
DiffEqBase.AbstractEnsembleSolution
EnsembleSolution
EnsembleTestSolution
EnsembleSummary
DiffEqBase.EnsembleAlgorithm
DiffEqBase.BasicEnsembleAlgorithm
DiffEqBase.AbstractEnsembleEstimator
EnsembleSerial
EnsembleDistributed
EnsembleThreads
EnsembleSplitThreads
```


## Functions

```@docs
DiffEqBase.calculate_ensemble_errors
DiffEqBase.batch_func
DiffEqBase.solve_batch
DiffEqBase.thread_monte
DiffEqBase.vector_batch_data_to_arr
```

### Analysis

```@docs
EnsembleAnalysis.get_timestep
EnsembleAnalysis.get_timepoint
EnsembleAnalysis.timestep_mean
EnsembleAnalysis.timestep_median
EnsembleAnalysis.timestep_meancor
EnsembleAnalysis.timestep_meancov
EnsembleAnalysis.timestep_meanvar
EnsembleAnalysis.timestep_quantile
EnsembleAnalysis.timestep_weighted_meancov
EnsembleAnalysis.timepoint_mean
EnsembleAnalysis.timepoint_median
EnsembleAnalysis.timepoint_meancor
EnsembleAnalysis.timepoint_meancov
EnsembleAnalysis.timepoint_meanvar
EnsembleAnalysis.timepoint_quantile
EnsembleAnalysis.timepoint_weighted_meancov
EnsembleAnalysis.timeseries_steps_mean
EnsembleAnalysis.timeseries_steps_median
EnsembleAnalysis.timeseries_steps_meancor
EnsembleAnalysis.timeseries_steps_meancov
EnsembleAnalysis.timeseries_steps_meanvar
EnsembleAnalysis.timeseries_steps_quantile
EnsembleAnalysis.timeseries_steps_weighted_meancov
EnsembleAnalysis.timeseries_point_mean
EnsembleAnalysis.timeseries_point_median
EnsembleAnalysis.timeseries_point_meancor
EnsembleAnalysis.timeseries_point_meancov
EnsembleAnalysis.timeseries_point_meanvar
EnsembleAnalysis.timeseries_point_quantile
EnsembleAnalysis.timeseries_point_weighted_meancov
EnsembleAnalysis.componentwise_mean
EnsembleAnalysis.componentwise_meancor
EnsembleAnalysis.componentwise_meancov
EnsembleAnalysis.componentwise_meanvar
EnsembleAnalysis.componentwise_vectors_timestep
EnsembleAnalysis.componentwise_weighted_meancov
EnsembleAnalysis.componentwise_vectors_timepoint
```


## Deprecated "Monte Carlo" aliases

The following are deprecated and are just aliases to the `Ensemble` equivalents:

* `AbstractMonteCarloProblem`: [`AbstractEnsembleProblem`](@ref DiffEqBase.AbstractEnsembleProblem)
* `MonteCarloProblem`: [`EnsembleProblem`](@ref)
* `MonteCarloAlgorithm`: [`EnsembleAlgorithm`](@ref DiffEqBase.EnsembleAlgorithm)
* `AbstractMonteCarloSolution`: [`AbstractEnsembleSolution`](@ref DiffEqBase.AbstractEnsembleSolution)
* `MonteCarloSolution`: [`EnsembleSolution`](@ref)
* `MonteCarloSummary`: [`EnsembleSummary`](@ref)
* `calculate_monte_errors`: [`calculate_ensemble_errors`](@ref DiffEqBase.calculate_ensemble_errors)
