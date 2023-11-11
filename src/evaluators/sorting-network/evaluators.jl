module SortingNetwork

export SortingOutcomeEvaluator

using Statistics
using CoEvo
using Random: AbstractRNG
using DataStructures: SortedDict
using CoEvo.Species: AbstractSpecies
using CoEvo.Individuals: Individual
using CoEvo.Evaluators: evaluate
using CoEvo.Evaluators.ScalarFitness: ScalarFitnessEvaluator, ScalarFitnessEvaluation, ScalarFitnessRecord
using ..Evaluators.Outcome: OutcomeScalarFitnessEvaluation

Base.@kwdef struct SortingOutcomeEvaluator <: CoEvo.Evaluators.Evaluator 
    max_sort_fitness::Float64
    maximize::Bool = true
    epsilon::Float64 = 1e-6
end

function CoEvo.Evaluators.evaluate(
    evaluator::SortingOutcomeEvaluator,
    ::AbstractRNG,
    species::AbstractSpecies,
    outcomes::Dict{Int, SortedDict{Int, Float64}}
) 
    # But we add a bonus for individuals that hit evaluator.max_sort_fitness
    # proportional to how close they are to the minimum number of codons
    # The bottom is the same as ScalarFitnessEvaluator
    individuals = [species.population ; species.children]
    filter!(individual -> individual.id in keys(outcomes), individuals)
    ids = [individual.id for individual in individuals]
    outcome_means = [mean(values(outcomes[id])) for id in ids]
    fitnesses = evaluator.maximize ? outcome_means : -outcome_means
    @assert evaluator.max_sort_fitness >= 0 "max_sort_fitness must be non-negative"
    min_fitness = minimum(fitnesses)
    @assert min_fitness >= 0 "Error: SortingOutcomeEvaluator encountered a negative fitness"
    # If individuals hit the maximum fitness, that means that they sorted all test cases correctly
    # We can then add the second fitness term, which is a bonus for having a small number of codons
    bonuses = [fit == evaluator.max_sort_fitness ? ind.genotype.max_codons - length(ind.genotype.codons) : 0
               for (ind, fit) in zip(individuals, fitnesses)]
    fitnesses = fitnesses + bonuses

    records = [ScalarFitnessRecord(id, fitness) for (id, fitness) in zip(ids, fitnesses)]
    sort!(records, by = x -> x.fitness, rev = true)

    evaluation = OutcomeScalarFitnessEvaluation(species.id, records, outcomes)
    return evaluation
end
end
