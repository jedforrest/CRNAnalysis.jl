module CRNAnalysis

# TODO consolidate
using RowEchelon, Catalyst, Graphs, Combinatorics, MetaGraphs, ModelingToolkit
using Catalyst, LinearAlgebra, InvertedIndices, SymPy, Symbolics
using Catalyst, Oscar, GAP, RowEchelon

include("polynomial_util.jl")
export symbolic_function

include("independent_decomposition.jl")
export independent_decomposition, construct_graph

include("generate_subnets.jl")
export generate_subnets, save_subnet, load_subnet

# TODO
include("rowspan_invariants.jl")
export rowspan_invariants

# TODO
include("get_oscar_steadystatesystem.jl")
export get_oscar_steadystatesystem, polynomial_ideal

end
