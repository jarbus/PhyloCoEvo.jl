using PhylogeneticTrees: PhylogeneticNode
using DataStructures: SortedDict
using CoEvo.Observers: Observation

struct QueueElement
    indA::PhylogeneticNode
    indB::PhylogeneticNode
    dist::Int
    caller::Union{QueueElement, Nothing}
    add_B::Bool
end


function find_k_nearest_interactions(
    ida::Int,
    idb::Int,
    psa::PhylogeneticSpecies,
    psb::PhylogeneticSpecies,
    individual_outcomes::Dict{Int, SortedDict{Int, Float64}},
    k::Int)
    """Find the k nearest interactions to `id1,id2` in `individual_outcomes`
    by searching over the trees in `ps1` and `ps2`
    """
    # (id1, id2, dist)
    k_nearest_interactions = Vector{Tuple{Int, Int, Int}}()

    n1 = psa.tree.tree[ida]
    n2 = psb.tree.tree[idb]

    # make empty queue element
    queue = [QueueElement(n1, n2, 0, nothing, true)]
    seen = Set{Tuple{Int, Int}}()
    iters = 0
    while length(queue) > 0
        el = popfirst!(queue)
        if el.dist == 10
            break
        end
        iters += 1
        seen_tup = (el.indA.id, el.indB.id)
        push!(seen, seen_tup)  # update seen set
        # println("popped ",el)
        if !isnothing(el.indA.parent)
            push!(queue, QueueElement(el.indA.parent, el.indB, el.dist+1, el, false))
        end
        if !isnothing(el.indB.parent)
            push!(queue, QueueElement(el.indA, el.indB.parent, el.dist+1, el, true))
        end
        for neiA in el.indA.children
            !isnothing(el.caller) && neiA == el.caller.indA && continue
            push!(queue, QueueElement(neiA, el.indB, el.dist+1, el, false))
        end
        el.add_B || continue
        for neiB in el.indB.children
            !isnothing(el.caller) &&  neiB == el.caller.indB && continue
            push!(queue, QueueElement(el.indA, neiB, el.dist+1, el, true))
        end
    end
    println("Iters: ", iters)
    println("Num_seen: ", length(seen))
    println("Num_dupes: ", iters - length(seen))
end


function estimate_outcomes!(
    individual_outcomes::Dict{Int, SortedDict{Int, Float64}},
    evaluators::Vector{<:Evaluator},
    species::Vector{<:AbstractSpecies},
    observations::Vector{<:Observation},)
    """Estimate outcomes for each individual in each species
    """
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
    println("len unevaluated_interactions: ", length(unevaluated_interactions))
    println("len all_interactions: ", length(all_interactions))
    println("len evaluated_interactions: ", length(evaluated_interactions))
end
