using DiffEqDocs, DifferentialEquations
using Test
using ExplicitImports

# Tests are just for docs generation

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(DiffEqDocs) === nothing
    @test check_no_stale_explicit_imports(DiffEqDocs) === nothing
end
