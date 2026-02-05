# TODO write tests

is_rational_function(f) = !isequal(denominator(f), 1)

implicit_form(eqn::Equation) = isequal(eqn.rhs, 0) ? eqn.lhs : (eqn.rhs - eqn.lhs)

# rational equation --> polynomial equation
function expand_rational_equation(eqn::Equation; simplify=true)
    f = implicit_form(eqn)
    # collect denominators
    denoms = Num[denominator(term)
        for term in Symbolics.terms(f)
        if is_rational_function(term)
    ]
    # expand and simplify
    new_f = f * prod(denoms)
    return simplify ? Symbolics.simplify(new_f) : new_f
end

function symbolic_function(eqn_rhs, vars, params)
    Symbolics.build_function(eqn_rhs, vars, params; expression=Val{false})[1]
end
