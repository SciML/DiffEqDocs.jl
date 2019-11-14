using Documenter,DiffEqBase,DiffEqProblemLibrary,DiffEqBiological

makedocs(modules=[DiffEqBase,DiffEqProblemLibrary,DiffEqBiological],
         doctest=false, clean=true,
         format = Documenter.HTML(analytics = "UA-90474609-3"),
         assets = ["assets/favicon.ico"],
         sitename="DifferentialEquations.jl",
         authors="Chris Rackauckas",
         pages = Any[
         "Home" => "index.md",
         "Tutorials" => Any[
           "tutorials/ode_example.md",
           "tutorials/advanced_ode_example.md",
           "tutorials/sde_example.md",
           "tutorials/rode_example.md",
           "tutorials/dde_example.md",
           "tutorials/dae_example.md",
           "tutorials/discrete_stochastic_example.md",
           "tutorials/jump_diffusion.md",
           "tutorials/bvp_example.md",
           "tutorials/additional.md"
         ],
         "Basics" => Any[
           "basics/overview.md",
           "basics/common_solver_opts.md",
           "basics/solution.md",
           "basics/plot.md",
           "basics/integrator.md",
           "basics/problem.md",
           "basics/faq.md",
           "basics/compatibility_chart.md"
         ],
         "Problem Types" => Any[
           "types/discrete_types.md",
           "types/ode_types.md",
           "types/dynamical_types.md",
           "types/split_ode_types.md",
           "types/steady_state_types.md",
           "types/bvp_types.md",
           "types/sde_types.md",
           "types/rode_types.md",
           "types/dde_types.md",
           "types/dae_types.md",
           "types/jump_types.md",
         ],
         "Solver Algorithms" => Any[
           "solvers/discrete_solve.md",
           "solvers/ode_solve.md",
           "solvers/dynamical_solve.md",
           "solvers/split_ode_solve.md",
           "solvers/steady_state_solve.md",
           "solvers/bvp_solve.md",
           "solvers/jump_solve.md",
           "solvers/sde_solve.md",
           "solvers/rode_solve.md",
           "solvers/dde_solve.md",
           "solvers/dae_solve.md",
           "solvers/benchmarks.md"
         ],
         "Additional Features" => Any[
           "features/performance_overloads.md",
           "features/diffeq_arrays.md",
           "features/diffeq_operator.md",
           "features/noise_process.md",
           "features/linear_nonlinear.md",
           "features/callback_functions.md",
           "features/callback_library.md",
           "features/ensemble.md",
           "features/io.md",
           "features/low_dep.md",
           "features/progress_bar.md"
         ],
         "Analysis Tools" => Any[
           "analysis/parameterized_functions.md",
           "analysis/parameter_estimation.md",
           "analysis/bifurcation.md",
           "analysis/sensitivity.md",
           "analysis/global_sensitivity.md",
           "analysis/uncertainty_quantification.md",
           "analysis/neural_networks.md",
           "analysis/dev_and_test.md"
         ],
         "Domain Modeling Tools" => Any[
             "models/multiscale.md",
             "models/physical.md",
             "models/financial.md",
             "models/biological.md",
             "models/external_modeling.md"
         ],
         "APIs" => Any[
             "DiffEqBase API" => [
                 "Overview" => "apis/diffeqbase/overview.md",
                 "apis/diffeqbase/functions.md",
                 "apis/diffeqbase/problems.md",
                 "apis/diffeqbase/solutions.md",
                 "apis/diffeqbase/solvers.md",
                 "apis/diffeqbase/de_types.md",
                 "apis/diffeqbase/operators.md",
                 "apis/diffeqbase/callbacks.md",
                 "apis/diffeqbase/interpolation.md",
                 "apis/diffeqbase/ensembles.md",
                 "apis/diffeqbase/data_arrays.md",
                 "apis/diffeqbase/noise.md",
                 "apis/diffeqbase/utility.md",
             ],
             "apis/diffeqbio.md"
         ],
         "Extra Details" => Any[
             "extras/timestepping.md",
         ]
         ])

deploydocs(
   repo = "github.com/JuliaDiffEq/DiffEqDocs.jl.git"
)
