# Integrators


## Types

```@docs
DiffEqBase.DEIntegrator
DiffEqBase.IntegratorIntervals
DiffEqBase.IntegratorTuples
```


## Interface

```@docs
DiffEqBase.initialize!(u, t, integrator::DiffEqBase.DEIntegrator, any_modified, c)
DiffEqBase.step!
DiffEqBase.addat!
DiffEqBase.get_tmp_cache
DiffEqBase.user_cache
DiffEqBase.u_cache
DiffEqBase.du_cache
DiffEqBase.ratenoise_cache
DiffEqBase.rand_cache
DiffEqBase.full_cache
DiffEqBase.resize_non_user_cache!
DiffEqBase.deleteat_non_user_cache!
DiffEqBase.addat_non_user_cache!
DiffEqBase.terminate!
DiffEqBase.get_du
DiffEqBase.get_du!
DiffEqBase.get_dt
DiffEqBase.get_proposed_dt
DiffEqBase.u_modified!
DiffEqBase.savevalues!
DiffEqBase.add_tstop!
DiffEqBase.add_saveat!
DiffEqBase.set_abstol!
DiffEqBase.set_reltol!
DiffEqBase.reinit!
DiffEqBase.auto_dt_reset!
DiffEqBase.change_t_via_interpolation!
DiffEqBase.addsteps!
DiffEqBase.reeval_internals_due_to_modification!
DiffEqBase.set_t!
DiffEqBase.set_u!
DiffEqBase.set_ut!
DiffEqBase.addat!
DiffEqBase.last_step_failed
DiffEqBase.check_error
DiffEqBase.check_error!
DiffEqBase.intervals
```


## Traits


```@docs
DiffEqBase.has_reinit
```
