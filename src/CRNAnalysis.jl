module CRNAnalysis

using RowEchelon, Catalyst, Graphs, Combinatorics, MetaGraphs, ModelingToolkit
# using Catalyst, Oscar, GAP, RowEchelon
# using Catalyst, LinearAlgebra, InvertedIndices, SymPy, Symbolics

include("independent_decomposition.jl")
export independent_decomposition, construct_graph

include("generate_subnets.jl")
export generate_subnets, save_subnet, load_subnet

include("get_Oscar_SteadyStatesystem.jl")

include("rowspan_polynomials.jl")

end
