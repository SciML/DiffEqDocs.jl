# Load OrdinaryDiffEq pages - must be available
ordinary_diffeq_pages_file = joinpath(@__DIR__, "ordinarydiffeq_pages.jl")
if !isfile(ordinary_diffeq_pages_file)
    error("OrdinaryDiffEq pages file not found at: $ordinary_diffeq_pages_file. Run the build process first.")
end

include(ordinary_diffeq_pages_file)

# Transform OrdinaryDiffEq pages to have the api/ordinarydiffeq prefix
function transform_ordinarydiffeq_pages(pages_array)
    transformed = []
    for page in pages_array
        if isa(page, String)
            push!(transformed, "api/ordinarydiffeq/" * page)
        elseif isa(page, Pair)
            key, value = page
            if isa(value, String)
                push!(transformed, key => "api/ordinarydiffeq/" * value)
            elseif isa(value, Vector)
                push!(transformed, key => transform_ordinarydiffeq_pages(value))
            end
        end
    end
    return transformed
end

ordinary_diffeq_pages = transform_ordinarydiffeq_pages(pages)

# Load StochasticDiffEq pages - if available
stochastic_diffeq_pages_file = joinpath(@__DIR__, "stochasticdiffeq_pages.jl")
stochastic_diffeq_pages = []
if isfile(stochastic_diffeq_pages_file)
    include(stochastic_diffeq_pages_file)

    # Transform StochasticDiffEq pages to have the api/stochasticdiffeq prefix
    function transform_stochasticdiffeq_pages(pages_array)
        transformed = []
        for page in pages_array
            if isa(page, String)
                push!(transformed, "api/stochasticdiffeq/" * page)
            elseif isa(page, Pair)
                key, value = page
                if isa(value, String)
                    push!(transformed, key => "api/stochasticdiffeq/" * value)
                elseif isa(value, Vector)
                    push!(transformed, key => transform_stochasticdiffeq_pages(value))
                end
            end
        end
        return transformed
    end

    stochastic_diffeq_pages = transform_stochasticdiffeq_pages(pages)
end

pages = Any["index.md",
    "getting_started.md",
    "Tutorials" => Any["tutorials/faster_ode_example.md",
        "tutorials/advanced_ode_example.md",
        "tutorials/sde_example.md",
        "tutorials/rode_example.md",
        "tutorials/dde_example.md",
        "tutorials/dae_example.md",
        "tutorials/jump_diffusion.md",
        "tutorials/bvp_example.md"],
    "Examples" => Any[
        "Beginner" => Any["examples/classical_physics.md",
            "examples/conditional_dosing.md",
            "examples/kepler_problem.md",
            "examples/outer_solar_system.md",
            "examples/min_and_max.md"],
        "Advanced" => Any["examples/spiking_neural_systems.md",
            "examples/beeler_reuter.md",
            "examples/diffusion_implicit_heat_equation.md"]],
    "Basics" => Any["basics/overview.md",
        "basics/common_solver_opts.md",
        "basics/solution.md",
        "basics/plot.md",
        "basics/integrator.md",
        "basics/problem.md",
        "basics/faq.md",
        "basics/compatibility_chart.md"],
    "Problem Types" => Any["types/discrete_types.md",
        "types/ode_types.md",
        "types/nonautonomous_linear_ode.md",
        "types/dynamical_types.md",
        "types/split_ode_types.md",
        "types/steady_state_types.md",
        "types/bvp_types.md",
        "types/sde_types.md",
        "types/sdae_types.md",
        "types/rode_types.md",
        "types/dde_types.md",
        "types/sdde_types.md",
        "types/dae_types.md"],
    "Solver Algorithms" => Any["solvers/discrete_solve.md",
        "solvers/ode_solve.md",
        "solvers/nonautonomous_linear_ode.md",
        "solvers/dynamical_solve.md",
        "solvers/split_ode_solve.md",
        "solvers/steady_state_solve.md",
        "solvers/bvp_solve.md",
        "solvers/sde_solve.md",
        "solvers/sdae_solve.md",
        "solvers/rode_solve.md",
        "solvers/dde_solve.md",
        "solvers/sdde_solve.md",
        "solvers/dae_solve.md",
        "solvers/benchmarks.md"],
    "Additional Features" => Any["features/performance_overloads.md",
        "features/diffeq_arrays.md",
        "features/diffeq_operator.md",
        "features/noise_process.md",
        "features/linear_nonlinear.md",
        "features/callback_functions.md",
        "features/callback_library.md",
        "features/ensemble.md",
        "features/io.md",
        "features/low_dep.md",
        "features/progress_bar.md"],
    "External Solver APIs" => Any["api/sundials.md",
        "api/daskr.md"],
    "OrdinaryDiffEq.jl API" => ordinary_diffeq_pages,
    "StochasticDiffEq.jl API" => stochastic_diffeq_pages,
    "Extra Details" => Any["extras/timestepping.md"]]
