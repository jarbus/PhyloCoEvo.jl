function find_k_nearest_interactions(
    ida::Int,
    idb::Int,
    pta::PhylogeneticTree,
    ptb::PhylogeneticTree,
    individual_outcomes::AbstractDict{Int, <:AbstractDict{Int, Float64}},
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
    @assert ida ∈ keys(pta.tree) "id1 $(ida) not in tree A"
    @assert idb ∈ keys(ptb.tree) "id2 $(idb) not in tree B"
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
    0 == length(k_nearest_interactions) && error("Found 0 interactions for $(ida),$(idb)")
    @assert length(seen) == iters "Seen set does not match number of iterations"
    return k_nearest_interactions
end

function compute_estimates(
    pairs::Vector{Tuple{Int, Int}},
    treeA::PhylogeneticTree,
    treeB::PhylogeneticTree,
    individual_outcomes::AbstractDict{Int, <:AbstractDict{Int, Float64}};
    k::Int,
    max_dist::Int)
    """For each pair of individuals in `pairs`, find the k nearest interactions
    in `individual_outcomes` and compute the weighted average outcome. 

    Returns:
    ========
    estimates: Vector{EstimatedOutcome}
        A vector of EstimatedOutcome objects for each pair of individuals in `pairs`
    """
    estimates = Vector{EstimatedOutcome}(undef, length(pairs))

    for (i, (ida, idb)) in enumerate(pairs)
        nearest = find_k_nearest_interactions(ida, idb, treeA, treeB, individual_outcomes, k, max_dist=max_dist)
        estimates[i] = EstimatedOutcome(ida, idb, nearest)
    end
    return estimates
end

function estimate!(
    estimator::PhylogeneticEstimator,
    individual_outcomes::Dict{Int, <:AbstractDict{Int, Float64}},
    speciesa::PhylogeneticSpecies,
    speciesb::PhylogeneticSpecies)
    @assert estimator.speciesa_id == speciesa.id 
    @assert estimator.speciesb_id == speciesb.id 

    # Merge all interactions into outcome cache
    two_layer_merge!(estimator.cached_outcomes, individual_outcomes)

    # Do we have sampled interactions?
    has_sampled_interactions =  speciesa.id ∈ keys(speciesb.randomly_sampled_interactions) &&
                                speciesb.id ∈ keys(speciesa.randomly_sampled_interactions) &&
                                (length(speciesa.randomly_sampled_interactions[speciesb.id]) > 0 ||
                                 length(speciesb.randomly_sampled_interactions[speciesa.id]) > 0)
    # If so, remove them from individual_outcomes and outcome cache
    if has_sampled_interactions
        sampled_interactions = union(
            [(i[1], i[2]) for i in speciesa.randomly_sampled_interactions[speciesb.id]],
            [(i[2], i[1]) for i in speciesb.randomly_sampled_interactions[speciesa.id]],
        )
        sampled_ids = Set(id for i in sampled_interactions for id in i)
        sampled_individual_outcomes = Dict(id=>SortedDict{Int,Float64}() for id in sampled_ids)
        for (id1, id2) in sampled_interactions
            @assert id1 in keys(individual_outcomes) "id1 $(id1) not in individual_outcomes"
            @assert id2 in keys(individual_outcomes[id1]) "id2 $(id2) not in individual_outcomes[$(id1)]"
            sampled_individual_outcomes[id1][id2] = individual_outcomes[id1][id2]
            delete!(individual_outcomes[id1], id2)
            delete!(estimator.cached_outcomes[id1], id2)

            sampled_individual_outcomes[id2][id1] = individual_outcomes[id2][id1]
            delete!(individual_outcomes[id2], id1)
            delete!(estimator.cached_outcomes[id2], id1)
        end
    end

    # Compute all unevaluated interactions between the two species
    evaluated_interactions = [
        (ind_i.id, ind_j.id)
        for ind_i in [speciesa.population; speciesa.children] 
            if ind_i.id ∈ keys(estimator.cached_outcomes)
        for ind_j in [speciesb.population; speciesb.children] 
            if ind_j.id ∈ keys(estimator.cached_outcomes) &&
               ind_j.id ∈ keys(estimator.cached_outcomes[ind_i.id])
    ]

    current_interactions = [
        (ind_i.id, ind_j.id)
        for ind_i in [speciesa.population; speciesa.children] 
            if ind_i.id ∈ keys(individual_outcomes)
        for ind_j in [speciesb.population; speciesb.children] 
            if ind_j.id ∈ keys(individual_outcomes[ind_i.id])
    ]

    all_interactions = [ 
        (ind_i.id, ind_j.id)
        for ind_i in [speciesa.population; speciesa.children]
        for ind_j in [speciesb.population; speciesb.children]]

    unevaluated_interactions = setdiff(all_interactions, evaluated_interactions)
    cached_interactions = setdiff(evaluated_interactions, current_interactions)

    # TODO add tests in test/estimator.jl
    
    # Compute estimates for sampled interactions
    if has_sampled_interactions
        # Estimate sampled interactions
        sample_estimates::Vector{EstimatedOutcome} = compute_estimates(
            sampled_interactions,
            speciesa.tree,
            speciesb.tree,
            estimator.cached_outcomes,
            k=estimator.k, max_dist=estimator.max_dist)

        sample_estimated_outcomes = estimates_to_outcomes(sample_estimates)
        two_layer_merge!(individual_outcomes, sample_estimated_outcomes, warn=true)

        # Compare sample estimates to actual outcomes
        esma, esmb = measure_estimation_samples(sample_estimates, sampled_individual_outcomes)

        # Add metrics to species
        if PhylogeneticEstimationSampleMeasurement ∉ keys(speciesa.measurements)
            speciesa.measurements[PhylogeneticEstimationSampleMeasurement] = Dict{String, Any}()
        end
        if PhylogeneticEstimationSampleMeasurement ∉ keys(speciesb.measurements)
            speciesb.measurements[PhylogeneticEstimationSampleMeasurement] = Dict{String, Any}()
        end
        speciesa.measurements[PhylogeneticEstimationSampleMeasurement][speciesb.id] = esma
        speciesb.measurements[PhylogeneticEstimationSampleMeasurement][speciesa.id] = esmb
    end

    # Compute estimates for the rest of the interactions
    estimates = compute_estimates(
                            unevaluated_interactions,
                            speciesa.tree,
                            speciesb.tree,
                            estimator.cached_outcomes,
                            k=estimator.k, max_dist=estimator.max_dist)
    estimated_individual_outcomes = estimates_to_outcomes(estimates)
    # merge estimated_individual_outcomes into individual_outcomes
    two_layer_merge!(individual_outcomes, estimated_individual_outcomes, warn=true)

    # Add cached interactions to individual_outcomes
    for (id1, id2) in cached_interactions
        individual_outcomes[id1][id2] = estimator.cached_outcomes[id1][id2]
        individual_outcomes[id2][id1] = estimator.cached_outcomes[id2][id1]
    end
end

function estimate!(estimator::PhylogeneticEstimator,
                   individual_outcomes::Dict{Int, <:AbstractDict{Int, Float64}},
                   species::Vector{<:AbstractSpecies})
    speciesa = find_species_by_id(estimator.speciesa_id, species)
    speciesb = find_species_by_id(estimator.speciesb_id, species)
    estimate!(estimator, individual_outcomes, speciesa, speciesb)
end
