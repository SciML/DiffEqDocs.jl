using Documenter,DiffEqDevTools,DiffEqBase,FiniteElementDiffEq,
      DiffEqProblemLibrary, StokesDiffEq

makedocs(modules=[DiffEqDevTools,DiffEqBase,FiniteElementDiffEq,
                  StokesDiffEq,OrdinaryDiffEq,DiffEqProblemLibrary],
         doctest=false, clean=true,
         format =:html,
         sitename="DifferentialEquations.jl",
         authors="Chris Rackauckas",
         pages = Any[
         "Home" => "index.md",
         "Tutorials" => Any[
           "tutorials/ode_example.md",
           "tutorials/sde_example.md",
           "tutorials/dae_example.md",
           "tutorials/fempoisson_example.md",
           "tutorials/femheat_example.md",
           "tutorials/femstochastic_example.md"
         ],
         "Basics" => Any[
           "basics/overview.md",
           "basics/common_solver_opts.md",
           "basics/solution.md",
           "basics/plot.md"
         ],
         "Equation Types" => Any[
           "types/ode_types.md",
           "types/sde_types.md",
           "types/fem_types.md",
           "types/stokes_types.md"
         ],
         "Solver Algorithms" => Any[
           "solvers/ode_solve.md",
           "solvers/sde_solve.md",
           "solvers/dae_solve.md",
           "solvers/fempoisson_solve.md",
           "solvers/femheat_solve.md",
           "solvers/fdmstokes_solve.md"
         ],
         "Additional Features" => Any[
           "features/performance_overloads.md",
           "features/callback_functions.md",
           "features/mesh.md",
           "features/output_specification.md",
           "features/conditional_dependencies.md",
           "features/progress_bar.md"
         ],
         "Analysis Tools" => Any[
           "analysis/parameterized_functions.md",
           "analysis/parameter_estimation.md",
           "analysis/sensitivity.md",
           "analysis/dev_and_test.md"
         ]
         ])

deploydocs(
   repo = "github.com/JuliaDiffEq/DiffEqDocs.jl.git",
   target = "build",
   osname = "linux",
   julia = "0.5",
   deps = nothing,
   make = nothing)
