module Lexicase

export LexicaseSelector

import CoEvo.Selectors: select

using Random: AbstractRNG
using DataStructures: OrderedDict
using CoEvo.Individuals: Individual
using PhyloCoEvo.Evaluators.Outcome: OutcomeScalarFitnessEvaluation
using CoEvo.Selectors: Selector

Base.@kwdef struct LexicaseSelector <: Selector
    n_parents::Int
end

function fast_max_filter!(source_idxs::Vector{Int},
        n_source_idxs::Int,
        target_idxs::Vector{Int},
        outcomes::Matrix{Float64},
        test_idx::Int)
    """Copies source_ids with the max outcome for a given test into the into first slots of target_ids."""
    cur_max = -Inf
    n_target_idxs = 0
    for i in 1:n_source_idxs
        @inbounds idx = source_idxs[i]
        @inbounds outcome = outcomes[idx, test_idx] # This is the lexicase selection bottleneck
        # Reset targets if we find new max
        if outcome > cur_max
            cur_max = outcome
            @inbounds target_idxs[1] = idx
            n_target_idxs = 1
        # Add to targets if we find another max
        elseif outcome == cur_max
            n_target_idxs += 1
            @inbounds target_idxs[n_target_idxs] = source_idxs[i]
        end
    end
    @assert cur_max != -Inf
    return n_target_idxs
end

function lexicase_select(
    rng::AbstractRNG,
    outcomes::Matrix{Float64})

    source_idxs = collect(1:size(outcomes, 1))
    target_idxs = Vector{Int}(undef, length(source_idxs))
    n_source_idxs = length(source_idxs)
    test_idxs = Set(1:size(outcomes, 2))

    while n_source_idxs > 1 && length(test_idxs) > 0
        # sample a test case without replacement
        rand_test_idx = pop!(test_idxs, rand(rng, test_idxs))
        # out of the remaining ids, choose the ones that have the best fitness on the test case
        n_source_idxs = fast_max_filter!(source_idxs, n_source_idxs, target_idxs, outcomes, rand_test_idx)
        # swap source and target ids
        source_idxs, target_idxs = target_idxs, source_idxs
    end
    # There is either one id left, or no test cases left. For both cases we can choose one at random
    @assert all(source_idxs[1:n_source_idxs] .<= size(outcomes, 1))
    rand(rng, source_idxs[1:n_source_idxs])
end

"""
    select(selector::LexicaseSelector, random_number_generator::AbstractRNG, new_population::Vector{<:Individual}, evaluation::OutcomeScalarFitnessEvaluation)
    
Perform lexicase selection to choose parents from a new population based on their fitness scores on multiple test cases.

# Arguments
- `selector::LexicaseSelector`: The lexicase selector object.
- `random_number_generator::AbstractRNG`: Random number generator.
- `new_population::Vector{<:Individual}`: The new population of individuals.
- `evaluation::OutcomeScalarFitnessEvaluation`: The outcome scalar fitness evaluation object.

# Returns
- `Vector{<:Individual}`: The selected parents.
"""

function select(
    selector::LexicaseSelector,
    random_number_generator::AbstractRNG, 
    new_population::Vector{<:I},
    evaluation::OutcomeScalarFitnessEvaluation
) where I<:Individual
    start = time()
    ids = collect(individual.id for individual in new_population) # index to id mapping
    test_cases = unique(test_id for id in ids for test_id in keys(evaluation.outcomes[id]))
    # initialize outcomes matrix to avoid repeated dict lookups
    outcomes = Matrix{Float64}(undef, length(ids), length(test_cases))
    for (idx, id) in enumerate(ids)
        for (test_idx, test_id) in enumerate(test_cases)
            outcomes[idx, test_idx] = evaluation.outcomes[id][test_id]
        end
    end
    parents = Vector{I}(undef, selector.n_parents)
    for i in 1:selector.n_parents
        new_parent_idx = lexicase_select(random_number_generator, outcomes)
        parents[i] = new_population[new_parent_idx]
    end
    endtime = time()
    println("Lexicase selection took $(endtime - start) seconds")
    return parents  
end

end
