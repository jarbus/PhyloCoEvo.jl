using PhylogeneticTrees: PhylogeneticNode
using DataStructures: SortedDict
using CoEvo.Observers: Observation


# TODO Metrics of interest
# Interaction distance DI between evaluated relatives and estimates (min, mean, max, std)
# Difference between estimates and evaluations for each unit of distance (min, mean, max, std)
# Total difference between estimates and evaluations (min, mean, max, std)
# Pairwise distances Dp between nodes in each population (min, mean, max, std)
# Number of interactions estimated
# Number of interactions evaluated


struct QueueElement
    indA::PhylogeneticNode
    indB::PhylogeneticNode
    dist::Int
    caller::Union{QueueElement, Nothing} # Prevent infinite loops by never returning to the caller
    add_B::Bool # Whether to add inds from tree B
end

struct RelatedOutcome
    """A related outcome is an interaction between two individuals, 
    and the distance to a query pair of individuals. Used for computing
    weighted averages of outcomes.
    """
    ida::Int
    idb::Int
    dist::Int
    outcomea::Float64
    outcomeb::Float64
end

function find_k_nearest_interactions(
    ida::Int,
    idb::Int,
    pta::PhylogeneticTree,
    ptb::PhylogeneticTree,
    individual_outcomes::Dict{Int, SortedDict{Int, Float64}},
    k::Int;
    max_dist::Int)
    """Find the k nearest interactions to `ida,idb` in `individual_outcomes`
    by searching over the trees in `pta` and `ptb`.

    This algorithm can also be thought of as a "Dual-Breadth-First-Search", in that
    it traverses two trees, one pair of nodes at a time. To avoid reaching the same
    pair from two different directions, we only search nodes in tree B when the
    current pair was added from a search in tree A.

    Arguments:
    =========
    id1, id2: Int
        The ids of the individuals to find the nearest interactions for
    pta, ptb: PhylogeneticTree
        The phylogenetic trees to search over. id1 must be in pta and id2 must be in ptb
    individual_outcomes: Dict{Int, SortedDict{Int, Float64}}
        A dictionary mapping individual ids to a sorted dictionary of interactions
        and their outcomes. Contains evaluated interactions.
    k: Int
        The number of nearest interactions to find

    Returns:
    ========
    k_nearest_interactions: Vector{RelatedOutcome}
        We return a vector of RelatedOutcome objects instead of a weighted average
        for testing purposes and code-reuse. We can compute different types of weighted
        averages from this vector.
    """
    k_nearest_interactions = Vector{RelatedOutcome}()

    # Initialize the search at the interaction
    n1 = pta.tree[ida]
    n2 = ptb.tree[idb]
    # Check if the interaction is already in the dictionary
    if  n2.id in keys(individual_outcomes[n1.id]) && 
        n1.id in keys(individual_outcomes[n2.id])
        error("Interaction $(n1.id),$(n2.id) already in dictionary")
    end
    queue = [QueueElement(n1, n2, 0, nothing, true)]
    seen = Set{Tuple{Int, Int}}()
    iters = 0
    # TODO: Profile and optimize
    while length(queue) > 0
        el = popfirst!(queue)
        el.dist > max_dist && break
        iters += 1
        ids = (el.indA.id, el.indB.id)
        push!(seen, ids)
        # If the interaction is in the set of outcomes, add it to the list and break if we have enough
        @assert ids[1] in keys(individual_outcomes) "id1 $(ids[1]) not in individual_outcomes"
        @assert ids[2] in keys(individual_outcomes) "id2 $(ids[2]) not in individual_outcomes"
        if ids[2] in keys(individual_outcomes[ids[1]]) && ids[1] in keys(individual_outcomes[ids[2]])
            outcomea = individual_outcomes[ids[1]][ids[2]]
            outcomeb = individual_outcomes[ids[2]][ids[1]]
            push!(k_nearest_interactions, RelatedOutcome(ids[1], ids[2], el.dist, outcomea, outcomeb))
            length(k_nearest_interactions) >= k && break
        end
        
        # Add parent interactions to queue first, as they are more likely to be found:
        # 1. Results from previous generations can be cached
        # 2. For "all vs best" or "all vs parents", we expect parent interactions to 
        #    occur more frequently than child interactions

        # Add Parent of a
        if (isnothing(el.caller) || el.indA.parent != el.caller.indA) && !isnothing(el.indA.parent)
            push!(queue, QueueElement(el.indA.parent, el.indB, el.dist+1, el, false))
        end
        # Add Parent of b if we are supposed to search on the B tree
        if (isnothing(el.caller) || el.indB.parent != el.caller.indB) && !isnothing(el.indB.parent) && el.add_B
            push!(queue, QueueElement(el.indA, el.indB.parent, el.dist+1, el, true))
        end

        # Add children of A
        for neiA in el.indA.children
            !isnothing(el.caller) && neiA == el.caller.indA && continue
            push!(queue, QueueElement(neiA, el.indB, el.dist+1, el, false))
        end
        # Add children of B if we are supposed to search on the B tree
        el.add_B || continue
        for neiB in el.indB.children
            !isnothing(el.caller) &&  neiB == el.caller.indB && continue
            push!(queue, QueueElement(el.indA, neiB, el.dist+1, el, true))
        end
    end
    if length(k_nearest_interactions) < k
        @warn "Found $(length(k_nearest_interactions)) < $k interactions for $(ida),$(idb)"
    end
    @assert length(seen) == iters "Seen set does not match number of iterations"
    return k_nearest_interactions
end

function weighted_average_outcome(related_outcomes::Vector{RelatedOutcome})
    """Compute the weighted average outcome of a set of related outcomes

    Arguments:
    =========
    related_outcomes: Vector{RelatedOutcome}

    Returns:
    ========
    weighted_average_a: Float64
        The weighted average outcome of the first individual
    weighted_average_b: Float64
        The weighted average outcome of the second individual
    """
    k = length(related_outcomes)
    if k == 1
        weights = [1.0]
    else
        dists = [related_outcomes[i].dist for i in 1:k]
        inv_dist =  sum(dists) .- dists
        weights = inv_dist ./ sum(inv_dist)
    end
    weighted_average_a, weighted_average_b = 0.0, 0.0
    for i in 1:k
        @inbounds weighted_average_a += related_outcomes[i].outcomea * weights[i]
        @inbounds weighted_average_b += related_outcomes[i].outcomeb * weights[i]
    end
    return weighted_average_a, weighted_average_b
end

function compute_estimates(
    pairs::Vector{Tuple{Int, Int}},
    treeA::PhylogeneticTree,
    treeB::PhylogeneticTree,
    individual_outcomes::Dict{Int, SortedDict{Int, Float64}};
    k::Int,
    max_dist::Int)
    """For each pair of individuals in `pairs`, find the k nearest interactions
    in `individual_outcomes` and compute the weighted average outcome. 

    TODO Save a fraction of estimates as samples to compute statistics of interest.

    TODO compute accuracy of estimates by comparing to evaluated interactions.

    Returns:
    ========
    estimates: Dict{Int, Dict{Int, Float64}}
        A dictionary mapping individuals to their weighted outcome scores against their opponents
    """
    # TODO: add ways to compute metrics of interest
    # get all unique ids in pairs
    ids = Set{Int}(p[1] for p in pairs) ∪ Set{Int}(p[2] for p in pairs)
    estimates = Dict{Int, Dict{Int, Float64}}(i=>Dict{Int, Float64}() for i in ids)

    for (ida, idb) in pairs
        nearest = find_k_nearest_interactions(ida, idb, treeA, treeB, individual_outcomes, k, max_dist=max_dist)
        weighted_average_a, weighted_average_b = weighted_average_outcome(nearest)
        estimates[ida][idb] = weighted_average_a
        estimates[idb][ida] = weighted_average_b
    end
    return estimates
end

function two_layer_merge!(d1::Dict{Int, <:AbstractDict}, d2::Dict{Int, <:AbstractDict})
    """Merge dictionary of dictionaries `d2` into `d1` by merging the inner dictionaries
    if the key is in both dictionaries, and adding the key+dict if it is not in `d1`."""
    # TODO: Profile and Optimize
    for id in keys(d2)
        # If id is in individual_outcomes, merge the two dictionaries using merge!
        if id ∈ keys(d1)
            merge!(d1[id], d2[id])
        else
            @warn "Estimating all outcomes for individual $id"
            d1[id] = d2[id]
        end
    end
end


function estimate_outcomes!(
    individual_outcomes::Dict{Int, SortedDict{Int, Float64}},
    species::Vector{<:AbstractSpecies},
    k::Int,
    max_dist::Int)

    """Fills individual outcomes for all pairs of individuals in `species` by
    estimating the outcomes of unevaluated interactions. Also computes metrics
    of interest for each species relating to these estimates.
    """
    # TODO: Profile and Optimize
    evaluated_interactions = [
        (ind_i, ind_j)
        for i in 1:(length(species)-1)
        for j in (i+1):length(species)
        for ind_i in [species[i].population; species[i].children] if ind_i.id in keys(individual_outcomes)
        for ind_j in [species[j].population; species[j].children] if ind_j.id in keys(individual_outcomes)
    ]
    all_interactions = [ 
        (ind_i, ind_j)
        for i in 1:(length(species)-1)
        for j in (i+1):length(species)
        for ind_i in [species[i].population; species[i].children]
        for ind_j in [species[j].population; species[j].children]]

    unevaluated_interactions = setdiff(all_interactions, evaluated_interactions)

    # compute estimates between each pair of species and merge into individual_outcomes
    for i in 1:(length(species)-1)
        for j in (i+1):length(species)
            estimated_individual_outcomes = compute_estimates(
                                    unevaluated_interactions,
                                    species[i].tree,
                                    species[j].tree,
                                    individual_outcomes,
                                    k=k, max_dist=max_dist)
            @info "len unevaluated_interactions:  $(length(unevaluated_interactions))"
            @info "len all_interactions:  $(length(all_interactions))"
            @info "len evaluated_interactions: $(length(evaluated_interactions))"
            # merge estimated_individual_outcomes into individual_outcomes
            two_layer_merge!(individual_outcomes, estimated_individual_outcomes)
        end
    end
end