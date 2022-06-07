using Documenter, DiffEqBase, SciMLBase, DiffEqProblemLibrary, OrdinaryDiffEq

ODEProblemLibrary = DiffEqProblemLibrary.ODEProblemLibrary
ODEProblemLibrary.importodeproblems()

SDEProblemLibrary = DiffEqProblemLibrary.SDEProblemLibrary
SDEProblemLibrary.importsdeproblems()

DDEProblemLibrary = DiffEqProblemLibrary.DDEProblemLibrary
DDEProblemLibrary.importddeproblems()

DAEProblemLibrary = DiffEqProblemLibrary.DAEProblemLibrary
DAEProblemLibrary.importdaeproblems()

include("pages.jl")

makedocs(modules=[DiffEqBase,SciMLBase,DiffEqProblemLibrary,ODEProblemLibrary,SDEProblemLibrary,DDEProblemLibrary,DAEProblemLibrary,OrdinaryDiffEq],
         doctest=false, clean=true,
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical="https://diffeq.sciml.ai/stable/"),
         sitename="DifferentialEquations.jl",
         authors="Chris Rackauckas",
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

deploydocs(
   repo = "github.com/SciML/DiffEqDocs.jl.git"
)
