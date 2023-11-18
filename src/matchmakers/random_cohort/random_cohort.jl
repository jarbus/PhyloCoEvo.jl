module RandomCohort
export RandomCohortMatchMaker, make_matches

using Random
using ...Utils: two_rand
using CoEvo
using CoEvo.MatchMakers: MatchMaker
using CoEvo.Matches.Basic: BasicMatch
using CoEvo.Species: AbstractSpecies
using CoEvo.MatchMakers.AllvsAll: AllvsAllMatchMaker, make_matches


Base.@kwdef mutable struct RandomCohortMatchMaker <: MatchMaker
    """Make cohorts of n_matches_per_ind individuals from each 
    species and perform all vs all in the cohorts. Does not support
    different n_matches_per_ind for 3+ species.

    TODO: support different n_matches_per_ind for 3+ species
    """
    n_matches_per_ind::Dict{String, Int}
    n_samples::Int
    cohorts::Vector{Symbol} = [:population, :children]
    gen::Int = 0
    n_gens_before_sampling::Int = 10
end


function get_individual_ids_from_cohorts(
    species::AbstractSpecies, matchmaker::RandomCohortMatchMaker
)
    individuals = vcat([getfield(species, cohort) for cohort in matchmaker.cohorts]...)
    ids = [individual.id for individual in individuals]
    return ids
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::RandomCohortMatchMaker, 
    rng::AbstractRNG,
    interaction_id::String, 
    species1::AbstractSpecies, 
    species2::AbstractSpecies
)
    """Create `matchmaker.n_matches` random matches between the two species. Log all interactions
    as randomly sampled in the phylogeneticspecies."""

    # We assume all pairwise matchups between species will happen each generation, so
    # we only increment the whenever we see a specific pair of species
    species_ids = [species1.id, species2.id]
    sorted_keys = sort(collect(keys(matchmaker.n_matches_per_ind)))
    if sorted_keys[1] ∈ species_ids && sorted_keys[2] ∈ species_ids
        matchmaker.gen += 1
    end
    n_samples = matchmaker.gen > matchmaker.n_gens_before_sampling ? matchmaker.n_samples : 0 

    if matchmaker.gen <= matchmaker.n_gens_before_sampling
        allvsall_mm = AllvsAllMatchMaker([:population,:children])
        return make_matches(allvsall_mm, rng, interaction_id, [species1, species2])
    end

    ids_1 = shuffle(rng, get_individual_ids_from_cohorts(species1, matchmaker))
    ids_2 = shuffle(rng, get_individual_ids_from_cohorts(species2, matchmaker))
    child_ids_1 = Set(ind.id for ind in species1.children)
    child_ids_2 = Set(ind.id for ind in species2.children)
    n_matches_per_ind1 = matchmaker.n_matches_per_ind[species1.id]
    n_matches_per_ind2 = matchmaker.n_matches_per_ind[species2.id]
    @assert length(ids_2) % n_matches_per_ind1 == 0 "species $(species2.id) of len $(length(ids_1)) not divisible by $(n_matches_per_ind1)"
    @assert length(ids_1) % n_matches_per_ind2 == 0 "species $(species1.id) of len $(length(ids_2)) not divisible by $(n_matches_per_ind2)"

    cohorts_1 = [ids_2[i:i+n_matches_per_ind1-1] 
                 for i in 1:n_matches_per_ind1:length(ids_2)]
    cohorts_2 = [ids_1[i:i+n_matches_per_ind2-1] 
                 for i in 1:n_matches_per_ind2:length(ids_1)]
    @assert length(cohorts_1) == length(cohorts_2) "cohorts of different length $(length(cohorts_1)) != $(length(cohorts_2))"
    @assert length(cohorts_1) > 1 || n_samples == 0 "n_samples > 0 but only one cohort"
    # shuffle cohorts
    cohorts_1 = shuffle!(rng, cohorts_1)
    cohorts_2 = shuffle!(rng, cohorts_2)
    match_ids = Set((id_1, id_2) for (ids_1, ids_2) in zip(cohorts_2, cohorts_1)
                                 for (id_1, id_2) in Iterators.product(ids_1, ids_2))
    if n_samples == 0
        matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
        return matches
    end

    # Add sampled interaction set
    if species2.id ∉ keys(species1.randomly_sampled_interactions)
        species1.randomly_sampled_interactions[species2.id] = Set{Tuple{Int,Int}}()
    end
    if species1.id ∉ keys(species2.randomly_sampled_interactions)
        species2.randomly_sampled_interactions[species1.id] = Set{Tuple{Int,Int}}()
    end

    children_matchups = Set((id_1, id_2) for (id_1, id_2) in vec(collect(Iterators.product(child_ids_1, child_ids_2))))
    non_matched_children_matchups = shuffle(rng, setdiff(children_matchups, match_ids) |> collect)
    @assert length(children_matchups) > 0 "No children matchups"
    @assert length(non_matched_children_matchups) > n_samples "length(non_matched_children_matchups) $(length(non_matched_children_matchups)) > n_samples $(n_samples)"

    # Add random samples to species if they exist
    non_cohort_samples1 = non_matched_children_matchups[1:n_samples]
    non_cohort_samples2 = [(id_2, id_1) for (id_1, id_2) in non_cohort_samples1]
    
    union!(species1.randomly_sampled_interactions[species2.id], non_cohort_samples1)
    union!(species2.randomly_sampled_interactions[species1.id], non_cohort_samples2)

    # Add random samples to match ids
    union!(match_ids, non_cohort_samples1)

    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
end


function CoEvo.MatchMakers.make_matches(
    matchmaker::RandomCohortMatchMaker,
    rng::AbstractRNG,
    interaction_id::String,
    all_species::Vector{<:AbstractSpecies}
)
    if length(all_species) != 2
        throw(ErrorException("Only two-entity interactions are supported for now."))
    end
    species1 = all_species[1]
    species2 = all_species[2]
    matches = CoEvo.MatchMakers.make_matches(matchmaker, rng, interaction_id, species1, species2)
    return matches
end

end
