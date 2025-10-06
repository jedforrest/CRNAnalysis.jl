# TODO approach to getting steady states depends on the deficiency of the subnetworks
# (see Example_Retina.ipynb)

function get_oscar_steadystatesystem(rn)
    #Set up the equations
    #INPUT: rn
    #OUTPUT: system of polynomial equations

    N = netstoichmat(rn) #stoichiometric matrix
    Y = substoichmat(rn) #reactant matrix, (columns correspond to exponent vectors of the monomials)
    W = Oscar.rref(matrix(QQ,conservationlaws(N)))[2]
    n, r = size(N) #number of species and reactions
    d = size(W)[1] #number of conservation laws


    # Set up the polynomial system
    T, k, c = polynomial_ring(QQ, "k"=>1:r,"c"=>1:d;internal_ordering=:deglex)
    K = fraction_field(T)
    B, x = polynomial_ring(K, "x"=>1:n;internal_ordering=:deglex)
    # polynomialSystem = Array{fmpq_mpoly, 1}()
    polynomialSystem = Vector{elem_type(B)}()

    #array for the equation system
    F = []

    #Create the monomial with exponent Y[:,j]
    monomials = []

    for j in 1:r
        mono = B(1)
        for i in 1:n
            mono = mono*x[i]^Y[i,j]
        end
        push!(monomials,mono)
    end

    #Create the steady state equations
    for i in 1:n
        f = B(0)
        for j in 1:r
            #Create the monomial with exponent Y[:,j]
            f = f + N[i,j]*k[j]*monomials[j]
        end
        push!(F,f)
    end


    return F
end

function polynomial_ideal(eqn_rhs, vars, params)
    # Convert Catalyst symbolic variables into Julia Symbol types
    varnames = tosymbol.(vars, escape=false)
    paramnames = tosymbol.(params)

    # Create polynomial ring in Oscar
    CC, oscar_coeffs = polynomial_ring(QQ, paramnames)
    ff = fraction_field(CC)
    RR, oscar_vars = polynomial_ring(ff, varnames)

    # Map Catalyst variables to Oscar variables
    cat_var_params = [vars; params]
    oscar_var_params = [oscar_vars; oscar_coeffs]
    cat_to_oscar = Dict(cat => oscar for (cat, oscar) in zip(cat_var_params, oscar_var_params))
    oscar_to_cat = Dict((oscar => cat) for (cat, oscar) in cat_to_oscar)

    # build Oscar polynomial by substituting oscar vars in Catalyst equations RHS (Right Hand Sides)
    polys = map(eqn_rhs) do rhs
        if rhs isa Number  # is a constant e.g. zero
            RR(rhs)
        else
            Symbolics.substitute(rhs, cat_to_oscar)
        end
    end
    ideal(RR, polys)
end


# # we may want to map back to Catalyst/Symbolics
# function to_symbolic_polynomial(poly::MPolyRingElem)
#     # create new symbolic vars from ring of poly
#     var_syms = Symbolics.variable.(gens(parent(poly)))
#     coeff_syms = Symbolics.variable.(gens(coefficient_ring(poly)))

#     # coefficients
#     coeff_terms = Oscar.coefficients(poly)
#     coeffs = zeros(Num, length(coeff_terms))
#     for (i, term) in enumerate(coeff_terms)
#         cf_and_es = coefficients_and_exponents(term.num)
#         coeff_polyn = sum(Int(cf) * prod(coeff_syms .^ es) for (cf, es) in cf_and_es)
#         coeffs[i] = coeff_polyn
#     end
#     # polynomial variables
#     exp_vecs = collect(Oscar.exponents(poly))
#     xs = [prod(var_syms .^ e_vec) for e_vec in exp_vecs]
#     coeffs' * xs
# end


function lift_symbolics(f, I, xs)
    A = coordinates(f, I)
    B = zeros(Num, length(A))
    for (i, aa) in enumerate(A)
        if length(aa.coeffs) == 0
            continue
        else
            ff_coeff = aa.coeffs[]
            cf_and_es = coefficients_and_exponents(ff_coeff.num)
            polyn = sum(Int(cf) * prod(xs .^ es) for (cf, es) in cf_and_es)
            B[i] = polyn
        end
    end
    return B
end
