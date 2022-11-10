using Documenter, DiffEqBase, SciMLBase, OrdinaryDiffEq
import ODEProblemLibrary, SDEProblemLibrary, DDEProblemLibrary, DAEProblemLibrary

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(modules = [
             DiffEqBase,
             SciMLBase,
             ODEProblemLibrary,
             SDEProblemLibrary,
             DDEProblemLibrary,
             DAEProblemLibrary,
             OrdinaryDiffEq,
         ],
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         doctest = false, clean = true,
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://diffeq.sciml.ai/stable/"),
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
