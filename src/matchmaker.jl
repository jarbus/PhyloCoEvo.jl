using Random
using CoEvo.MatchMakers: MatchMaker
using CoEvo.Matches.Basic: BasicMatch

Base.@kwdef struct PhylogeneticMatchMaker <: MatchMaker
    cohorts::Vector{Symbol} = [:population, :children]
end

Base.@kwdef struct RandomCohortMatchMaker <: MatchMaker
    n_matches_per_ind::Dict{String, Int}
    cohorts::Vector{Symbol} = [:population, :children]
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::PhylogeneticMatchMaker, 
    rng::AbstractRNG,
    interaction_id::String, 
    species1::PhylogeneticSpecies, 
    species2::PhylogeneticSpecies
)
    # TODO: update phylogenetic species
    ids_1 = get_individual_ids_from_cohorts(species1, matchmaker)
    ids_2 = get_individual_ids_from_cohorts(species2, matchmaker)

    # Make sure that there are no previous interactions in the species
    @assert length(species1.randomly_sampled_interactions) == 0
    @assert length(species2.randomly_sampled_interactions) == 0
    match_ids = vec(collect(Iterators.product(ids_1, ids_2)))
    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
end

function get_individual_ids_from_cohorts(
    species::AbstractSpecies, matchmaker::Union{RandomCohortMatchMaker, PhylogeneticMatchMaker}
)
    individuals = vcat([getfield(species, cohort) for cohort in matchmaker.cohorts]...)
    ids = [individual.id for individual in individuals]
    return ids
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::RandomCohortMatchMaker, 
    rng::AbstractRNG,
    interaction_id::String, 
    species1::PhylogeneticSpecies, 
    species2::PhylogeneticSpecies
)
    """Create `matchmaker.n_matches` random matches between the two species. Log all interactions
    as randomly sampled in the phylogeneticspecies."""
    ids_1 = get_individual_ids_from_cohorts(species1, matchmaker)
    ids_2 = get_individual_ids_from_cohorts(species2, matchmaker)
    n_matches_per_ind1 = matchmaker.n_matches_per_ind[species1.id]
    n_matches_per_ind2 = matchmaker.n_matches_per_ind[species2.id]
    @assert length(ids_2) % n_matches_per_ind1 == 0 \
        "species $(species2.id) of len $(length(ids_1)) not divisible by $(n_matches_per_ind1)"
    @assert length(ids_1) % n_matches_per_ind2 == 0 \
        "species $(species1.id) of len $(length(ids_2)) not divisible by $(n_matches_per_ind2)"

    cohorts_1 = [ids_2[i:i+n_matches_per_ind1-1] 
                 for i in 1:n_matches_per_ind1:length(ids_2)]
    cohorts_2 = [ids_1[i:i+n_matches_per_ind2-1] 
                 for i in 1:n_matches_per_ind2:length(ids_1)]
    # shuffle cohorts
    cohorts_1 = shuffle!(rng, cohorts_1)
    cohorts_2 = shuffle!(rng, cohorts_2)
    # create matches between
    # create groups of ids_1 and ids_2 of size n_matches_per_ind
    # Make sure that there are no previous interactions in the species
    @assert length(species1.randomly_sampled_interactions) == 0
    @assert length(species2.randomly_sampled_interactions) == 0
    species1.randomly_sampled_interactions = Set(match_ids)
    species2.randomly_sampled_interactions = Set(match_ids)
    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
end


function CoEvo.MatchMakers.make_matches(
    matchmaker::Union{PhylogeneticMatchMaker, RandomCohortMatchMaker},
    rng::AbstractRNG,
    interaction_id::String,
    all_species::Vector{PhylogeneticSpecies},
)
    if length(all_species) != 2
        throw(ErrorException("Only two-entity interactions are supported for now."))
    end
    species1 = all_species[1]
    species2 = all_species[2]
    matches = make_matches(matchmaker, rng, interaction_id, species1, species2)
    return matches
end
