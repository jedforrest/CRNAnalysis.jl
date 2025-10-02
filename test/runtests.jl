using CRNAnalysis
using Catalyst
using Graphs
using Test

@testset "CRNAnalysis.jl" begin
    @testset "Decomposition" include("decomposition_tests.jl")

end
