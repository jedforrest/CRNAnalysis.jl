function generate_subnets(reaction_network::ReactionSystem, subnet_indices; verbose=true)

    subnets = ReactionSystem[]

    for (i, idxs) in enumerate(subnet_indices)
        verbose && println("Subnetwork $i: ", idxs)

        name = "$(reaction_network.name)_SN$i"
        subnet = ReactionSystem(reactions(reaction_network)[idxs], default_t(); name=Symbol(name))

        push!(subnets, subnet)
    end
    subnets
end
function generate_subnets(reaction_network::ReactionSystem; verbose=true)
    subnet_indices = independent_decomposition(reaction_network; verbose)
    generate_subnets(reaction_network, subnet_indices; verbose)
end

function save_subnet(rn::ReactionSystem; path=".")
    mkpath(path)
    filepath = joinpath(path, "$(rn.name).jl")
    save_reactionsystem(filepath, rn)
    filepath
end

function load_subnet(net_name; path=".")
    include("$path/$net_name.jl")
end
