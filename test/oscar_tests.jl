@testset "Symbolic conversion" begin
    #TODO
end

@testset "Oscar steady state equations" begin
    f1 = get_oscar_steadystatesystem(retina_sn1)
    I1 = ideal(f1)

    x = gens(base_ring(I1))
    k = gens(coefficient_ring(I1))
    @test I1[1] == -k[3]*x[1]*x[3] + k[1]
    @test I1[2] == 0
    @test I1[3] == -k[3]*x[1]*x[3] + k[2]*x[2]

    G1 = Oscar.groebner_basis(I1, ordering=lex([x[1],x[3],x[2]]),  complete_reduction = true)
    @test G1[1] == k[2]*x[2] - k[1]
    @test G1[2] == k[3]*x[1]*x[3] - k[1]
end

@testset "Groebner basis" begin
    #TODO

    # setpoint
end


# I = polynomial_ideal(eqn_rhs, vars, params)

# to_symbolic_polynomial.(gens(G))

# params = parameters(rn)
# sym_coeffs = lift_symbolics(G[1], I, params)
