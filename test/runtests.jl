using CRNAnalysis
using Test
using Catalyst
using Graphs
using Oscar
using Symbolics

retina_rn = include("../data/retina.jl")
retina_sn1 = include("../data/retina_sn1.jl")

@testset "CRNAnalysis.jl" begin
    @testset "Independent Decomposition" include("decomposition_tests.jl")
    @testset "Rowspan Invariants" include("invariant_tests.jl")
    @testset "Oscar" include("oscar_tests.jl")
end
