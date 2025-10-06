test_equal(expr1, expr2) = isequal(simplify(expr1 - expr2), 0)
test_equal_vec(vec1, vec2) = length(vec1) != length(vec2) ? false : all(test_equal.(vec1, vec2))

@testset "Rowspan polynomials" begin
    @unpack θ, μ, φ, λ, β, κ, Ce, Ce1, Ce2 = retina_rn
    sn1, sn2, sn3 = generate_subnets(retina_rn)[1:3]; # TODO don't need to load all subnets
    complex_ids = [1, 3]

    b, invariant = rowspan_invariants(sn1, complex_ids)
    rpa_eqn = invariant * massactionvector(sn1)
    @test test_equal_vec(b, [-1/θ, 0, 1/θ])
    @test test_equal_vec(invariant, [-μ/θ 0 1 0 0])
    @test isequal(rpa_eqn[], -μ/θ + Ce)

    _, invariant = rowspan_invariants(sn2, complex_ids)
    rpa_eqn = invariant * massactionvector(sn2)
    @test test_equal_vec(invariant, [-φ/λ 0 1 0 0])
    @test isequal(rpa_eqn[], -φ/λ + Ce1)

    _, invariant = rowspan_invariants(sn3, complex_ids)
    rpa_eqn = invariant * massactionvector(sn3)
    @test test_equal_vec(invariant, [-β/κ 0 1 0 0])
    @test isequal(rpa_eqn[], -β/κ + Ce2)
end

# TODO all subnetworks
sn = generate_subnets(retina_rn)
complex_ids = [1, 3]
b, invariant = rowspan_invariants(sn[4], complex_ids)
b
invariant
