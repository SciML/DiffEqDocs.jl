using Documenter, DiffEqBase, SciMLBase, OrdinaryDiffEq, OrdinaryDiffEqCore
import ODEProblemLibrary,
       SDEProblemLibrary, DDEProblemLibrary, DAEProblemLibrary, BVProblemLibrary
using Sundials, DASKR

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

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
        OrdinaryDiffEqCore,
        Sundials, DASKR
    ],
    linkcheck = true,
    linkcheck_ignore = ["https://www.izhikevich.org/publications/spikes.htm",
        "https://biojulia.net/post/hardware/",
        "https://archimede.dm.uniba.it/~testset/report/pollu.pdf",
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
        "https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html"
    ],
    doctest = false, clean = true,
    warnonly = [:missing_docs],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DiffEqDocs/stable/"),
    sitename = "DifferentialEquations.jl",
    authors = "Chris Rackauckas",
    pages = pages)

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
