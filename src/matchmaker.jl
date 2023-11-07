using CoEvo.MatchMakers: MatchMaker
using CoEvo.Matches.Basic: BasicMatch
using CoEvo.MatchMakers.AllvsAll: get_individual_ids_from_cohorts

Base.@kwdef struct PhylogeneticMatchMaker <: MatchMaker
    cohorts::Vector{Symbol} = [:population, :children]
end

Base.@kwdef struct RandomMatchMaker <: MatchMaker
    n_matches::Int
    cohorts::Vector{Symbol} = [:population, :children]
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::PhylogeneticMatchMaker, 
    interaction_id::String, 
    species1::PhylogeneticSpecies, 
    species2::PhylogeneticSpecies
)
    # TODO: update phylogenetic species
    ids_1 = get_individual_ids_from_cohorts(species1, matchmaker)
    ids_2 = get_individual_ids_from_cohorts(species2, matchmaker)
    match_ids = vec(collect(Iterators.product(ids_1, ids_2)))
    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::PhylogeneticMatchMaker,
    rng::AbstractRNG,
    interaction_id::String,
    all_species::Vector{PhylogeneticSpecies},
)
    if length(all_species) != 2
        throw(ErrorException("Only two-entity interactions are supported for now."))
    end
    species1 = all_species[1]
    species2 = all_species[2]
    matches = make_matches(matchmaker, interaction_id, species1, species2)
    return matches
end
