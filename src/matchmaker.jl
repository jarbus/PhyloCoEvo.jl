export ParentsVsChildrenMatchMaker, RandomCohortMatchMaker, make_matches
using Random
using CoEvo.MatchMakers: MatchMaker
using CoEvo.Matches.Basic: BasicMatch

Base.@kwdef struct ParentsVsChildrenMatchMaker <: MatchMaker
    """Make matches between parents and children. All individuals
    in species.population are considered parents, and all individuals
    in species.children are considered children.
    
    Arguments:
        n_random_samples::Int: Number of additional random child v child 
        matches to make, used for estimating the accuracy of estimates.
    """
    n_samples::Int
end

Base.@kwdef struct RandomCohortMatchMaker <: MatchMaker
    """Make cohorts of n_matches_per_ind individuals from each 
    species and perform all vs all in the cohorts. Does not support
    different n_matches_per_ind for 3+ species.

    TODO: support different n_matches_per_ind for 3+ species
    """
    n_matches_per_ind::Dict{String, Int}
    n_samples::Int
    cohorts::Vector{Symbol} = [:population, :children]
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::ParentsVsChildrenMatchMaker, 
    rng::AbstractRNG,
    interaction_id::String, 
    species1::PhylogeneticSpecies, 
    species2::PhylogeneticSpecies
)
    parent1_ids = [ind.id for ind in species1.population]
    parent2_ids = [ind.id for ind in species2.population]
    child1_ids = [ind.id for ind in species1.children]
    child2_ids = [ind.id for ind in species2.children]
    
    parents_v_children_matches = Iterators.flatten(
                (Iterators.product(parent1_ids, child2_ids), 
                 Iterators.product(child1_ids, parent2_ids))
               ) |> collect |> vec

    # choose n_random_samples random child v child matches
    child_v_child_matches = Iterators.product(child1_ids, child2_ids) |> collect |> vec
    # shuffle matches
    shuffle!(rng, child_v_child_matches)
    # Choose n_random_samples random matches
    random_child_v_child_matches = Set(child_v_child_matches[1:matchmaker.n_samples])
    @assert length(random_child_v_child_matches) == matchmaker.n_samples
    # Update phylogenetic species
    union!(species1.randomly_sampled_interactions, random_child_v_child_matches)
    union!(species2.randomly_sampled_interactions, random_child_v_child_matches)

    # Combine all matches
    match_ids = union(Set(parents_v_children_matches), Set(random_child_v_child_matches))
    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
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
    species1::PhylogeneticSpecies, 
    species2::PhylogeneticSpecies
)
    """Create `matchmaker.n_matches` random matches between the two species. Log all interactions
    as randomly sampled in the phylogeneticspecies."""
    ids_1 = get_individual_ids_from_cohorts(species1, matchmaker)
    ids_2 = get_individual_ids_from_cohorts(species2, matchmaker)
    n_matches_per_ind1 = matchmaker.n_matches_per_ind[species1.id]
    n_matches_per_ind2 = matchmaker.n_matches_per_ind[species2.id]
    @assert length(ids_2) % n_matches_per_ind1 == 0 "species $(species2.id) of len $(length(ids_1)) not divisible by $(n_matches_per_ind1)"
    @assert length(ids_1) % n_matches_per_ind2 == 0 "species $(species1.id) of len $(length(ids_2)) not divisible by $(n_matches_per_ind2)"

    cohorts_1 = [ids_2[i:i+n_matches_per_ind1-1] 
                 for i in 1:n_matches_per_ind1:length(ids_2)]
    cohorts_2 = [ids_1[i:i+n_matches_per_ind2-1] 
                 for i in 1:n_matches_per_ind2:length(ids_1)]
    @assert length(cohorts_1) == length(cohorts_2) "cohorts of different length $(length(cohorts_1)) != $(length(cohorts_2))"
    @assert length(cohorts_1) > 1 || matchmaker.n_samples == 0 "n_samples > 0 but only one cohort"
    # shuffle cohorts
    cohorts_1 = shuffle!(rng, cohorts_1)
    cohorts_2 = shuffle!(rng, cohorts_2)
    match_ids = Set((id_1, id_2) for (ids_1, ids_2) in zip(cohorts_1, cohorts_2)
                    for (id_1, id_2) in vec(collect(Iterators.product(ids_1, ids_2))))
    
    # Choose n_random_samples between members of non-corresponding cohorts
    non_cohort_samples = Set()
    for i in 1:matchmaker.n_samples
        c1, c2 = two_rand(rng, 1:length(cohorts_1))
        id_1 = rand(rng, cohorts_1[c1])
        id_2 = rand(rng, cohorts_2[c2])
        push!(non_cohort_samples, (id_1, id_2))
    end

    # Update phylogenetic species
    union!(species1.randomly_sampled_interactions, non_cohort_samples)
    union!(species2.randomly_sampled_interactions, non_cohort_samples)

    # Add random samples to match ids
    union!(match_ids, non_cohort_samples)

    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
end


function CoEvo.MatchMakers.make_matches(
    matchmaker::Union{ParentsVsChildrenMatchMaker, RandomCohortMatchMaker},
    rng::AbstractRNG,
    interaction_id::String,
    all_species::Vector{<:PhylogeneticSpecies{<:Individual}}
)
    if length(all_species) != 2
        throw(ErrorException("Only two-entity interactions are supported for now."))
    end
    species1 = all_species[1]
    species2 = all_species[2]
    matches = make_matches(matchmaker, rng, interaction_id, species1, species2)
    return matches
end
