using DiffEqDocs, DifferentialEquations
using Test
using ExplicitImports

# Tests are just for docs generation

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqDocs) === nothing
    @test check_no_stale_explicit_imports(DiffEqDocs) === nothing
end

@testset "Documenter errors" begin
    makefile = readlines(joinpath(pkgdir(DiffEqDocs), "docs", "make.jl"))
    warnonly = only(filter(line -> occursin("warnonly =", line), makefile))
    @test !occursin(":cross_references", warnonly)

    imports = join(makefile, '\n')
    @test occursin("OrdinaryDiffEqFunctionMap", imports)
    @test occursin("OrdinaryDiffEqNewmark", imports)
end
