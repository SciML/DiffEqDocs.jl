using Documenter,DiffEqBase,DiffEqPDEBase,DiffEqProblemLibrary

makedocs(modules=[DiffEqBase,DiffEqPDEBase,DiffEqProblemLibrary],
         doctest=false, clean=true,
         format =:html,
         sitename="DifferentialEquations.jl",
         authors="Chris Rackauckas",
         pages = Any[
         "Home" => "index.md",
         "Tutorials" => Any[
           "tutorials/ode_example.md",
           "tutorials/sde_example.md",
           "tutorials/rode_example.md",
           "tutorials/dde_example.md",
           "tutorials/dae_example.md",
           "tutorials/discrete_stochastic_example.md",
           "tutorials/jump_diffusion.md",
           "tutorials/fempoisson_example.md",
           "tutorials/femheat_example.md",
           "tutorials/femstochastic_example.md"
         ],
         "Basics" => Any[
           "basics/overview.md",
           "basics/common_solver_opts.md",
           "basics/solution.md",
           "basics/plot.md",
           "basics/integrator.md",
           "basics/compatibility_chart.md"
         ],
         "Equation Types" => Any[
           "types/discrete_types.md",
           "types/ode_types.md",
           "types/refined_ode_types.md",
           "types/sde_types.md",
           "types/rode_types.md",
           "types/dde_types.md",
           "types/dae_types.md",
           "types/refined_dae_types.md",
           "types/jump_types.md",
           "types/fem_types.md",
         ],
         "Solver Algorithms" => Any[
           "solvers/discrete_solve.md",
           "solvers/ode_solve.md",
           "solvers/refined_ode_solve.md",
           "solvers/sde_solve.md",
           "solvers/rode_solve.md",
           "solvers/dde_solve.md",
           "solvers/dae_solve.md",
           "solvers/fempoisson_solve.md",
           "solvers/femheat_solve.md"
         ],
         "Additional Features" => Any[
           "features/performance_overloads.md",
           "features/diffeq_arrays.md",
           "features/noise_process.md",
           "features/linear_nonlinear.md",
           "features/callback_functions.md",
           "features/callback_library.md",
           "features/monte_carlo.md",
           "features/low_dep.md",
           "features/mesh.md",
           "features/progress_bar.md"
         ],
         "Analysis Tools" => Any[
           "analysis/parameterized_functions.md",
           "analysis/parameter_estimation.md",
           "analysis/bifurcation.md",
           "analysis/sensitivity.md",
           "analysis/uncertainty_quantification.md",
           "analysis/dev_and_test.md"
         ],
         "Modeling Tools" => Any[
             "models/multiscale.md",
             "models/financial.md",
             "models/biological.md"
         ]
         ])

deploydocs(
   repo = "github.com/JuliaDiffEq/DiffEqDocs.jl.git",
   target = "build",
   osname = "linux",
   julia = "0.5",
   deps = nothing,
   make = nothing)
