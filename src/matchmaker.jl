using CoEvo.Ecosystems.Interactions.MatchMakers.Abstract: MatchMaker
# using CoEvo.Ecosystems.Species.Abstract: AbstractSpecies
# using Random: AbstractRNG

using CoEvo.Ecosystems.Interactions.MatchMakers.Matches.Basic: BasicMatch
import CoEvo.Ecosystems.Interactions.MatchMakers.Interfaces: make_matches

"""
    PhylogeneticMatchMaker <: MatchMaker

Struct for specifying the matchmaking strategy between two species. 
Supports two types:
- `:plus`: Matches all combinations of individuals, including both current population and children.
- `:comma`: Matches only current population or children, depending on which is non-empty.

# Fields
- `type::Symbol`: Specifies the type of matchmaking to be done (`:plus` or `:comma`).
"""
Base.@kwdef struct PhylogeneticMatchMaker <: MatchMaker 
    type::Symbol = :plus
end

"""
    (mm::PhylogeneticMatchMaker)(sp1::Species, sp2::Species) -> Vector{Tuple{Int, Int}}

Matchmaking function to create pairs of individuals between two species based on the type of `PhylogeneticMatchMaker`.

# Arguments
- `mm::PhylogeneticMatchMaker`: Matchmaking configuration.
- `sp1::Species`: First species for matchmaking.
- `sp2::Species`: Second species for matchmaking.

# Returns
- A vector of pairs (Tuples), each pair represents a match between an individual from `sp1` and `sp2`.

# Errors
- Throws an error if an invalid type is set in `PhylogeneticMatchMaker`.
"""
function make_matches(
    matchmaker::PhylogeneticMatchMaker, 
    interaction_id::String, 
    sp1::AbstractSpecies, 
    sp2::AbstractSpecies
)
    if matchmaker.type == :comma
        ids1 = length(sp1.children) == 0 ? collect(keys(sp1.pop)) : collect(keys(sp1.children))
        ids2 = length(sp2.children) == 0 ? collect(keys(sp2.pop)) : collect(keys(sp2.children))
    elseif matchmaker.type == :plus
        ids1 = [collect(keys(sp1.pop)); collect(keys(sp1.children))]
        ids2 = [collect(keys(sp2.pop)); collect(keys(sp2.children))]
    else
        error("Invalid PhylogeneticMatchMaker type: $(matchmaker.type)")
    end
    match_ids = vec(collect(Iterators.product(ids1, ids2)))
    matches = [
        BasicMatch(interaction_id, [id1, id2]) for (id1, id2) in match_ids
    ]
    return matches
end

function make_matches(
    matchmaker::PhylogeneticMatchMaker,
    ::AbstractRNG,
    interaction_id::String,
    all_species::Dict{String, <:AbstractSpecies},
    species_ids::Vector{String}
)
    if length(species_ids) != 2
        throw(ErrorException("Only two-entity interactions are supported for now."))
    end
    species1 = all_species[species_ids[1]]
    species2 = all_species[species_ids[2]]
    matches = make_matches(matchmaker, interaction_id, species1, species2)
    return matches
end