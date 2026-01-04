using Documenter, DiffEqBase, SciMLBase, OrdinaryDiffEq, OrdinaryDiffEqBDF,
    OrdinaryDiffEqCore, StochasticDiffEq, DelayDiffEq, SteadyStateDiffEq, DiffEqCallbacks
import ODEProblemLibrary,
    SDEProblemLibrary, DDEProblemLibrary, DAEProblemLibrary, BVProblemLibrary
using Sundials, DASKR

# Use development versions for API documentation
import Pkg
Pkg.develop("OrdinaryDiffEq")
Pkg.develop("StochasticDiffEq")

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

# Copy OrdinaryDiffEq.jl documentation
ordinartdiffeq_docs_root = joinpath(dirname(pathof(OrdinaryDiffEq)), "..", "docs")
ordinartdiffeq_docs_path = joinpath(ordinartdiffeq_docs_root, "src")
if isdir(ordinartdiffeq_docs_path)
    # Create the OrdinaryDiffEq API directory in the docs
    ordinary_diffeq_dest = joinpath(@__DIR__, "src", "api", "ordinarydiffeq")
    mkpath(dirname(ordinary_diffeq_dest))

    # Copy all the docs from OrdinaryDiffEq.jl
    cp(ordinartdiffeq_docs_path, ordinary_diffeq_dest, force = true)

    # Copy the pages.jl file from OrdinaryDiffEq.jl
    ordinary_diffeq_pages_dest = joinpath(@__DIR__, "ordinarydiffeq_pages.jl")
    ordinary_diffeq_pages_file = joinpath(ordinartdiffeq_docs_root, "pages.jl")
    cp(ordinary_diffeq_pages_file, ordinary_diffeq_pages_dest, force = true)

    # Copy the common_first_steps.jl and common_imex_first_steps.jl files from OrdinaryDiffEq.jl
    for fname in ["common_first_steps.jl", "common_imex_first_steps.jl"]
        common_first_steps_dest = joinpath(@__DIR__, fname)
        common_first_steps_file = joinpath(ordinartdiffeq_docs_root, fname)
        cp(common_first_steps_file, common_first_steps_dest, force = true)
    end
end

# Copy StochasticDiffEq.jl documentation
stochasticdiffeq_docs_root = joinpath(dirname(pathof(StochasticDiffEq)), "..", "docs")
stochasticdiffeq_docs_path = joinpath(stochasticdiffeq_docs_root, "src")
if isdir(stochasticdiffeq_docs_path)
    # Create the StochasticDiffEq API directory in the docs
    stochastic_diffeq_dest = joinpath(@__DIR__, "src", "api", "stochasticdiffeq")
    mkpath(dirname(stochastic_diffeq_dest))

    # Copy all the docs from StochasticDiffEq.jl
    cp(stochasticdiffeq_docs_path, stochastic_diffeq_dest, force = true)

    # Copy the pages.jl file from StochasticDiffEq.jl
    stochastic_diffeq_pages_dest = joinpath(@__DIR__, "stochasticdiffeq_pages.jl")
    stochastic_diffeq_pages_file = joinpath(stochasticdiffeq_docs_root, "pages.jl")
    cp(stochastic_diffeq_pages_file, stochastic_diffeq_pages_dest, force = true)
end

ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

include("pages.jl")

makedocs(
    modules = [
        DiffEqBase,
        SciMLBase,
        ODEProblemLibrary,
        SDEProblemLibrary,
        DDEProblemLibrary,
        DAEProblemLibrary,
        BVProblemLibrary,
        OrdinaryDiffEq,
        OrdinaryDiffEq.OrdinaryDiffEqAdamsBashforthMoulton,
        OrdinaryDiffEq.OrdinaryDiffEqBDF,
        OrdinaryDiffEq.OrdinaryDiffEqDefault,
        OrdinaryDiffEq.OrdinaryDiffEqExplicitRK,
        OrdinaryDiffEq.OrdinaryDiffEqExponentialRK,
        OrdinaryDiffEq.OrdinaryDiffEqExtrapolation,
        OrdinaryDiffEq.OrdinaryDiffEqFeagin,
        OrdinaryDiffEq.OrdinaryDiffEqFIRK,
        OrdinaryDiffEq.OrdinaryDiffEqHighOrderRK,
        OrdinaryDiffEq.OrdinaryDiffEqIMEXMultistep,
        OrdinaryDiffEq.OrdinaryDiffEqLinear,
        OrdinaryDiffEq.OrdinaryDiffEqLowOrderRK,
        OrdinaryDiffEq.OrdinaryDiffEqLowStorageRK,
        OrdinaryDiffEq.OrdinaryDiffEqNordsieck,
        OrdinaryDiffEq.OrdinaryDiffEqPDIRK,
        OrdinaryDiffEq.OrdinaryDiffEqPRK,
        OrdinaryDiffEq.OrdinaryDiffEqQPRK,
        OrdinaryDiffEq.OrdinaryDiffEqRKN,
        OrdinaryDiffEq.OrdinaryDiffEqRosenbrock,
        OrdinaryDiffEq.OrdinaryDiffEqSDIRK,
        OrdinaryDiffEq.OrdinaryDiffEqSSPRK,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedRK,
        OrdinaryDiffEq.OrdinaryDiffEqSymplecticRK,
        OrdinaryDiffEq.OrdinaryDiffEqTsit5,
        OrdinaryDiffEq.OrdinaryDiffEqVerner,
        StochasticDiffEq,
        DelayDiffEq,
        SteadyStateDiffEq,
        DiffEqCallbacks,
        Sundials, DASKR,
    ],
    remotes = nothing,
    linkcheck = true,
    linkcheck_ignore = [
        "https://www.izhikevich.org/publications/spikes.htm",
        "https://biojulia.net/post/hardware/",
        "https://archimede.dm.uniba.it/~testset/report/pollu.pdf",
        r"https://archimede.uniba.it/~bvpsolvers/testsetbvpsolvers/\?page_id=\d+",
        "http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Demos_Pitagora/DemoHires/demohires.pdf",
        "https://www.radford.edu/%7Ethompson/RP/nonnegative.pdf",
        "http://www.radford.edu/~thompson/vodef90web/problems/demosnodislin/Demos_Pitagora/DemoOrego/demoorego.pdf",
        "https://zenodo.org/record/5883749#.Yg-d698xmu4",
        "https://www.sciencedirect.com/science/article/abs/pii/S0375960109009591",
        "https://www.sciencedirect.com/science/article/abs/pii/0375960176901018",
        "https://www.worldscientific.com/doi/abs/10.1142/S0218127499001024",
        "https://www.worldscientific.com/doi/abs/10.1142/S0218127499001383",
        "https://www.wolframalpha.com/input/?i=u%27%3D-sqrt%28u%29",
        "https://www.mathworks.com/help/simulink/gui/absolutetolerance.html",
        "https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html",
        "https://journals.ametsoc.org/view/journals/atsc/20/2/1520-0469_1963_020_0130_dnf_2_0_co_2.xml",
        "https://mathematica.stackexchange.com/questions/40122/help-to-plot-poincar%C3%A9-section-for-double-pendulum",
        "https://github.com/SciML/DiffEqProblemLibrary.jl/blob/master/lib/SDEProblemLibrary/src/SDEProblemLibrary.jl",
        "https://github.com/SciML/ColPrac/blob/master/README.md",
        "https://github.com/SciML/DiffEqDevTools.jl/blob/master/src/ode_tableaus.jl",
        "https://github.com/SciML/DiffEqProblemLibrary.jl/blob/master/lib/BVProblemLibrary/src/BVProblemLibrary.jl",
        "https://github.com/SciML/DiffEqProblemLibrary.jl/blob/master/lib/DDEProblemLibrary/src/DDEProblemLibrary.jl",
        "https://docs.sciml.ai/DiffEqDocs/stable/features/dae_initialization/",
    ],
    doctest = false, clean = true,
    warnonly = [:missing_docs, :docs_block],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DiffEqDocs/stable/",
        size_threshold = 500 * 2^10
    ),
    sitename = "DifferentialEquations.jl",
    authors = "Chris Rackauckas",
    pages = pages
)

#Redirect old links
# cd(joinpath(@__DIR__, "build")) do
#     for (root, dirs, files) in walkdir(".")
#         for file in files
#             path = relpath(joinpath(root, file), ".")
#             m = match(r"(.+)/index\.html$", path)
#             m === nothing && continue
#             redirect = "$(m[1]).html"
#             @info "Adding redirect for $(m[1]) from $(redirect)"
#             isfile(redirect) && (@warn "$redirect exists, skip"; continue)
#             open(redirect, "w") do io
#                 write(io, """
#                 <meta http-equiv="refresh" content="0; url=$(basename(m[1]))/"/>
#                 """)
#             end
#         end
#     end
# end

deploydocs(repo = "github.com/SciML/DiffEqDocs.jl.git")
