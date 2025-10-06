let
# Serialised using Catalyst version v1.11.7.

# Independent variable:
@parameters t

# Parameters:
ps = @parameters φ λ ψ

# Species:
sps = @species Sci1(t) Ce1(t) C1(t)

# Reactions:
rxs = [
	Reaction(φ, nothing, [Sci1], nothing, [1]),
	Reaction(λ, [Ce1], [Ce1, C1], [1], [1, 1]),
	Reaction(ψ, [Sci1, C1], nothing, [1, 1], nothing)
]

# Declares ReactionSystem model:
ReactionSystem(rxs, t, sps, ps; name = :Retina_SN2)

end