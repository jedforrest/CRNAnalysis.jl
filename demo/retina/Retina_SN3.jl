let
# Serialised using Catalyst version v1.11.7.

# Independent variable:
@parameters t

# Parameters:
ps = @parameters β κ ε

# Species:
sps = @species Sci2(t) Ce2(t) C2(t)

# Reactions:
rxs = [
	Reaction(β, nothing, [Sci2], nothing, [1]),
	Reaction(κ, [Ce2], [Ce2, C2], [1], [1, 1]),
	Reaction(ε, [Sci2, C2], nothing, [1, 1], nothing)
]

# Declares ReactionSystem model:
ReactionSystem(rxs, t, sps, ps; name = :Retina_SN3)

end