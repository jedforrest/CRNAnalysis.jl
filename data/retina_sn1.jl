let
# Serialised using Catalyst version v1.11.7.

# Independent variable:
@parameters t

# Parameters:
ps = @parameters μ θ η

# Species:
sps = @species Sci(t) Ce(t) C(t)

# Reactions:
rxs = [
	Reaction(μ, nothing, [Sci], nothing, [1]),
	Reaction(θ, [Ce], [Ce, C], [1], [1, 1]),
	Reaction(η, [Sci, C], nothing, [1, 1], nothing)
]

# Declares ReactionSystem model:
ReactionSystem(rxs, t, sps, ps; name = :Retina_SN1)

end