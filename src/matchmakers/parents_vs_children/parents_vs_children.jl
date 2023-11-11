module ParentsVsChildren
export ParentsVsChildrenMatchMaker, make_matches

using Random
using CoEvo
using CoEvo.MatchMakers: MatchMaker
using CoEvo.MatchMakers.AllvsAll: AllvsAllMatchMaker, make_matches

using CoEvo.Matches.Basic: BasicMatch
using CoEvo.Species: AbstractSpecies

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


function CoEvo.MatchMakers.make_matches(
    matchmaker::ParentsVsChildrenMatchMaker, 
    rng::AbstractRNG,
    interaction_id::String, 
    species1::AbstractSpecies, 
    species2::AbstractSpecies
)
    # If children have no parents, then run all vs all with no samples
    c_parents = [p for s in [species1, species2] for ind in s.children for p in ind.parent_ids]
    n_expected_parents = length(species1.children) + length(species2.children)
    @assert length(c_parents) == 0 || length(c_parents) == n_expected_parents "We expect either 0 parents or $n_expected_parents parents, but got $(length(c_parents)) parents."
    if length(c_parents) == 0
        allvsall_mm = AllvsAllMatchMaker([:population,:children])
        return make_matches(allvsall_mm, rng, interaction_id, [species1, species2])
    end
    
    # Otherwise, match parents v children
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

function CoEvo.MatchMakers.make_matches(
    matchmaker::ParentsVsChildrenMatchMaker,
    rng::AbstractRNG,
    interaction_id::String,
    all_species::Vector{<:AbstractSpecies}
)
    if length(all_species) != 2
        throw(ErrorException("Only two-entity interactions are supported for now."))
    end
    species1 = all_species[1]
    species2 = all_species[2]
    matches = make_matches(matchmaker, rng, interaction_id, species1, species2)
    return matches
end

end
