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

"""
    lexicase_select(rng::AbstractRNG, all_ids::Set{Int}, all_test_cases::Set{Int}, outcomes::AbstractDict{Int,AbstractDict{Int,Float64}})
    
Perform lexicase selection to choose an individual from a population based on the outcomes on multiple test cases.

# Arguments
- `rng::AbstractRNG`: Random number generator.
- `all_ids::Set{Int}`: Set of all individual IDs in the population.
- `all_test_cases::Set{Int}`: Set of all test case IDs.
- `outcomes::AbstractDict{Int,AbstractDict{Int,Float64}}`: Outcome scores of all individuals on all test cases.

# Returns
- `Int`: ID of the selected individual.
"""
function lexicase_select(
    rng::AbstractRNG,
    all_ids::Set{Int},
    all_test_cases::Set{Int},
    outcomes::AbstractDict{Int,<:AbstractDict{Int,Float64}})
    ids = copy(all_ids)
    test_cases = copy(all_test_cases)
    while length(ids) > 1 && length(test_cases) > 0
        # sample a test case without replacement
        rand_test_case = pop!(test_cases, rand(rng, test_cases))
        # out of the remaining ids, choose the ones that have the best fitness on the test case
        best_outcome = maximum(outcomes[id][rand_test_case] for id in ids)
        ids = filter(id -> outcomes[id][rand_test_case] == best_outcome, ids)
    end
    # There is either one id left, or no test cases left. For both cases we can choose one at random
    rand(rng, ids)
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
    
    ids = Set(individual.id for individual in new_population)
    id_to_index = Dict(new_population[i].id => i for i in eachindex(new_population))
    test_cases = Set(test_id for id in ids for test_id in keys(evaluation.outcomes[id]))
    parents = I[]
    for i in 1:selector.n_parents
        new_parent_id = lexicase_select(random_number_generator, ids, test_cases, evaluation.outcomes)
        push!(parents, new_population[id_to_index[new_parent_id]])
    end
    return parents  
end

end
