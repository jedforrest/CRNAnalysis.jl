let
# Serialised using Catalyst version v1.11.7.

# Independent variable:
@parameters t

# Parameters:
ps = @parameters p11 k5 p5 α ω k9 p12 k12 k22 k24 k6 p6 k21 p15

# Species:
sps = @species Cp1(t) CO(t) Cp(t) Ce(t) Ce1(t) H(t) H1(t) HR(t) HR1(t) Cg(t) B(t)

# Reactions:
rxs = [
	Reaction(p11, [Cp1], [CO], [1], [1]),
	Reaction(k5, [Cp], [Ce], [1], [1]),
	Reaction(p5, [Cp1], [Ce1], [1], [1]),
	Reaction(α, nothing, [H], nothing, [1]),
	Reaction(ω, nothing, [H1], nothing, [1]),
	Reaction(k9, [H, HR], [HR, Ce], [1, 1], [1, 1]),
	Reaction(p12, [H1, HR1], [HR1, Ce1], [1, 1], [1, 1]),
	Reaction(k12, [Cp, CO], [Cg], [1, 1], [1]),
	Reaction(k22, [Ce, Cg], [B], [1, 1], [1]),
	Reaction(k24, [B], nothing, [1], nothing),
	Reaction(k6, [Ce], [Cp], [1], [1]),
	Reaction(p6, [Ce1], [Cp1], [1], [1]),
	Reaction(k21, [Ce], nothing, [1], nothing),
	Reaction(p15, [Ce1], nothing, [1], nothing)
]

# Declares ReactionSystem model:
ReactionSystem(rxs, t, sps, ps; name = :Retina_SN23)

end