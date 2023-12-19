module RandomParentsVsChildren
export RandomParentsVsChildrenMatchMaker, make_matches

using Random
using CoEvo
using CoEvo.MatchMakers: MatchMaker
using CoEvo.MatchMakers.AllvsAll: AllvsAllMatchMaker, make_matches

using CoEvo.Matches.Basic: BasicMatch
using ...Species.Phylogenetic: PhylogeneticSpecies

mutable struct UniqueRandomSelector
    elements::Vector
    shuffled_indices::Vector
end

UniqueRandomSelector(elements::Vector) =
    UniqueRandomSelector(elements, shuffle(collect(1:length(elements))))

function get_random_element(selector::UniqueRandomSelector)
    if isempty(selector.shuffled_indices)
        selector.shuffled_indices = shuffle(collect(1:length(selector.elements)))
    end
    return selector.elements[pop!(selector.shuffled_indices)]
end

Base.@kwdef struct RandomParentsVsChildrenMatchMaker <: MatchMaker
    """Make matches between parents and children, sometimes making child v child matches. All individuals in species.population are
    considered parents, and all individuals in species.children are considered
    children. We expect all children to have approximately the same number of
    matches, and all parents to have approximately the same number of matches.
    
    Arguments:
        n_samples::Int: Number of additional random child v child 
        matches to make, used for estimating the accuracy of estimates.
        n_cvc_per_child::Int: Number of parent v child matchups to instead run as child v child matchups.
    """
    n_samples::Int
    n_cvc_per_child::Int
end


function CoEvo.MatchMakers.make_matches(
    matchmaker::RandomParentsVsChildrenMatchMaker, 
    rng::AbstractRNG,
    interaction_id::String, 
    species1::PhylogeneticSpecies, 
    species2::PhylogeneticSpecies
)
    # If children have no parents, then run all vs all with no samples
    c_parents = [p for s in [species1, species2] for ind in s.children for p in ind.parent_ids]
    n_expected_parents = length(species1.children) + length(species2.children)
    @assert length(c_parents) == 0 || 
        length(c_parents) == n_expected_parents "We expect either 0 parents or $n_expected_parents parents, but got $(length(c_parents)) parents."
    @assert length(c_parents) == 0 || 
        (matchmaker.n_cvc_per_child < length(species1.population) &&
         matchmaker.n_cvc_per_child < length(species2.population)) "We expect n_cvc_per_child to be less than the number of parents."


        
    if length(c_parents) == 0
        allvsall_mm = AllvsAllMatchMaker([:population,:children])
        return make_matches(allvsall_mm, rng, interaction_id, [species1, species2])
    end
    
    # Otherwise, match parents v children
    parent1_ids = [ind.id for ind in species1.population]
    parent2_ids = [ind.id for ind in species2.population]
    child1_ids = [ind.id for ind in species1.children]
    child2_ids = [ind.id for ind in species2.children]
    # print all ids
    parent1_selector = UniqueRandomSelector(parent1_ids)
    parent2_selector = UniqueRandomSelector(parent2_ids)
    child2_selector = UniqueRandomSelector(child2_ids)

    # Keep a set of all child matches, so we can make sure to not duplicate them when sampling
    child_v_child_matches = Iterators.product(child1_ids, child2_ids) |> collect |> Set

    match_pool = Set{Tuple{Int,Int}}()
    # Go over all children and match them with parents and children
    # First, do species 1
    for child1_id in child1_ids
        # Match children with parents of species 2
        for _ in 1:(length(parent2_ids) - matchmaker.n_cvc_per_child)
            parent2_id = get_random_element(parent2_selector)
            while (child1_id, parent2_id) in match_pool
                parent2_id = get_random_element(parent2_selector)
            end
            push!(match_pool, (child1_id, parent2_id))
            delete!(child_v_child_matches, (child1_id, parent2_id)) # Don't sample this match again
        end
        # Match with children with children of species 2
        for _ in 1:matchmaker.n_cvc_per_child
            child2_id = get_random_element(child2_selector)
            while (child1_id, child2_id) in match_pool
                child2_id = get_random_element(child2_selector)
            end
            push!(match_pool, (child1_id, child2_id))
            delete!(child_v_child_matches, (child1_id, child2_id)) # Don't sample this match again
        end
    end
    # Next, add parents of species 1 v children of species 2. We don't need to
    # repeat the child v child matches.
    for child2_id in child2_ids
        for _ in 1:(length(parent1_ids) - matchmaker.n_cvc_per_child)
            parent1_id = get_random_element(parent1_selector)
            while (parent1_id, child2_id) in match_pool
                parent1_id = get_random_element(parent1_selector)
            end
            push!(match_pool, (parent1_id, child2_id))
            delete!(child_v_child_matches, (parent1_id, child2_id)) # Don't sample this match again
        end
    end

    # from the remaining child v child matches, choose n_samples matches to use for estimating accuracy
    remaining_child_v_child_matches = collect(child_v_child_matches)
    shuffle!(rng, remaining_child_v_child_matches)
    remaining_child_v_child_matches = remaining_child_v_child_matches[1:matchmaker.n_samples]
    # Choose n_random_samples random matches
    random_child_v_child_matches1 = Set(remaining_child_v_child_matches)
    random_child_v_child_matches2 = Set((m[2],m[1]) for m in random_child_v_child_matches1)
    @assert length(random_child_v_child_matches1) == matchmaker.n_samples

    if species1 isa PhylogeneticSpecies && species2 isa PhylogeneticSpecies
        # Update phylogenetic species
        if species2.id ∉ keys(species1.randomly_sampled_interactions)
            species1.randomly_sampled_interactions[species2.id] = Set{Tuple{Int,Int}}()
        end
        if species1.id ∉ keys(species2.randomly_sampled_interactions)
            species2.randomly_sampled_interactions[species1.id] = Set{Tuple{Int,Int}}()
        end
        union!(species1.randomly_sampled_interactions[species2.id], random_child_v_child_matches1)
        union!(species2.randomly_sampled_interactions[species1.id], random_child_v_child_matches2)
    else
        @assert matchmaker.n_samples == 0 "We expect n_samples to be 0 for non-phylogenetic species."
    end
    # Combine all matches
    match_ids = union(
        match_pool,
        Set(random_child_v_child_matches1),
    )
    matches = [BasicMatch(interaction_id, [id_1, id_2]) for (id_1, id_2) in match_ids]
    return matches
end

function CoEvo.MatchMakers.make_matches(
    matchmaker::RandomParentsVsChildrenMatchMaker,
    rng::AbstractRNG,
    interaction_id::String,
    all_species::Vector{<:PhylogeneticSpecies}
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
