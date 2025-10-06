function independent_decomposition(reaction_network::ReactionSystem; verbose=true)

    # LINEAR COMBINATIONS
    netstoich = netstoichmat(reaction_network)
    linear_combo, pivots = linear_combinations(netstoich)

    ### GRAPH CONSTRUCTION
    G = construct_graph(linear_combo, pivots)  # FIXME rename

    # determine whether a finer decomposition exists, if not, end here
    if Graphs.is_connected(G)
        @info "No finer decomposition exists"
        return reaction_network
    end

    # determine the number of connected components of G
    # this is the number of partitions R will be decomposed to
    components = Graphs.connected_components(G)
    verbose && println("Number of independent subnetworks: $(length(components))")

    # PARTITION
    # initialize the vector of partitions
    P = []

    # basis vectors: assign them first into their respective partition
    for comp in components
        push!(P, pivots[comp])
    end

    # nonbasis vectors: they go to the same partition as the basis vectors
    # that form their linear combination
    for i = eachindex(P)
        for j = eachindex(P[i])

            # get the column number representing the basis vectors in 'linear_combo'
            col = findall(x -> x == P[i][j], pivots)[1]

            # check which reactions used a particular basis vector and assign them
            # to their respective partition
            append!(P[i], findall(x -> x != 0, linear_combo[:, col]))
        end
        # get only unique elements in each partition
        P[i] = unique(P[i])
    end

    parts = P

    # check that all reactions have been assigned a partition
    num_rxns = size(netstoich, 2)
    num_parts = length(reduce(vcat, parts))
    if num_parts != num_rxns
        missing_rxns = num_rxns - num_parts
        @warn "$missing_rxns reactions have been omitted in the partition"
        # TODO add this step
        # should repeat from the linear combination step (need test example)
    end

    parts  # TODO should return partition RN here
end

function construct_graph(linear_combo, pivots)

    # initialise an undirected graph G
    G = MetaGraph()  # TODO does this need to be a metagraph?

    # add vertices to G which are the vectors in basis
    for i = eachindex(pivots)
        Graphs.add_vertex!(G)
        name = string("R", pivots[i])
        set_prop!(G, i, :name, name)
    end

    abs_linear_combo = abs.(linear_combo)
    get_edges = first.(Tuple.(findall(x -> x > 1, sum(abs_linear_combo, dims=2))))

    #  initialize an array for sets of vertices that will form the edges
    vertex_sets = Vector[]

    # identify which vertices form edges in each reaction:
    # get those with non-zero coefficients in the linear combinations
    for edge in get_edges
        push!(vertex_sets, first.(Tuple.(findall(!iszero, linear_combo[edge, :]))))
    end

    # get all possible combinations (not permutations) of the reactions
    # involved in the linear combinations
    edges = collect.(Combinatorics.combinations.(vertex_sets, 2))

    #  get just the unique edges (if edges is non-empty)
    if !isempty(edges)
        edges = unique(reduce(vcat, edges))
    end

    # add these edges to graph G
    for edge in edges
        Graphs.add_edge!(G, edge[1], edge[2])
    end

    G
end
function construct_graph(reaction_network::ReactionSystem)
    netstoich = netstoichmat(reaction_network)
    linear_combo, pivots = linear_combinations(netstoich)
    construct_graph(linear_combo, pivots)
end

function linear_combinations(netstoich)

    num_rxns = size(netstoich, 2)

    # get the transpose of the stoichiometric matrix
    R = transpose(netstoich)

    # compute a maximal linearly independent set of vectors in R
    _, pivots = rref_with_pivots(netstoich)

    basis = R[pivots,:]

    # initialise matrix of linear combinations (#rxns x #basis elements)
    linear_combo = zeros(num_rxns, length(pivots))

    # get matrix which shows the linear combinations of basis rxns which
    # give the nonbasis rxns (basis rxns will have a row of zeros)
    for i = 1:num_rxns
        if !in(i, pivots)
            linear_combo[i,:] = transpose(basis) \ R[i,:]
        end
    end

    linear_combo = round.(linear_combo, digits=1)

    linear_combo, pivots
end
