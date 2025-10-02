retina_rn = include("../data/retina.jl")

# Test:
# - independent_decomposition
# - ind return early (connected graph)


@testset "construct_graph" begin

    netstoich = netstoichmat(retina_rn)
    linear_combo, pivots = CRNAnalysis.linear_combinations(netstoich)

    @test size(linear_combo) == (70, 45)
    @test length(pivots) == 45
    @test all(iszero(linear_combo[pivot_row, :]) for pivot_row in pivots)

    G = construct_graph(retina_rn)

    @test Graphs.nv(G) == 45
    @test Graphs.ne(G) == 50
end


@testset "independent_decomposition" begin

    ind_decomp = independent_decomposition(retina_rn; verbose=false)

    @test length(ind_decomp) == 27
    @test ind_decomp[1] == [1, 4, 7]
    @test maximum(length.(ind_decomp)) == 14

    # test that all reactions are correctly partitioned
    @test length(reduce(vcat, ind_decomp)) == numreactions(retina_rn)
end


@testset "generate_subnets" begin

    subnets = generate_subnets(retina_rn, verbose=false)
    sn1 = subnets[1]

    @test length(subnets) == 27
    @test sn1.name == :Retina_SN1
    @test numspecies(sn1) == 3
    @test deficiency(sn1) == 1

    # test that no further decomposition is needed
    out = independent_decomposition(sn1)
    @test out == sn1
end

@testset "save/load subnet files" begin

    temp_folder = mktempdir()

    subnets = generate_subnets(retina_rn, verbose=false)

    save_subnet.(subnets; path=temp_folder)

    # load the first and last subnets
    Retina_SN1 = load_subnet("Retina_SN1"; path=temp_folder)
    Retina_SN27 = load_subnet("Retina_SN27"; path=temp_folder)

    @test Retina_SN1 == subnets[1]
    @test Retina_SN27 == subnets[end]
end
