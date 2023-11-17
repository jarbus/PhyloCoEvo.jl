module SortingNetwork

export SortingNetworkEvaluator

using Statistics
using CoEvo
using Random: AbstractRNG
using DataStructures: SortedDict
using CoEvo.Species: AbstractSpecies
using CoEvo.Individuals: Individual
using CoEvo.Evaluators: evaluate, Evaluator, Evaluation
using CoEvo.Evaluators.ScalarFitness: ScalarFitnessEvaluator, ScalarFitnessEvaluation, ScalarFitnessRecord
using ..Evaluators.Outcome: OutcomeScalarFitnessEvaluation
using ...Genotypes.SortingNetwork: SortingNetworkGenotype

Base.@kwdef struct SortingNetworkEvaluator <: CoEvo.Evaluators.Evaluator 
    evaluator::Evaluator
    num_tests_per_parasite::Float64
end

function add_bonus_fitness(
    evaluator::SortingNetworkEvaluator,
    records::Vector{R},
    bonuses::Vector{Float64},
    outcomes::AbstractDict{Int, <:AbstractDict{Int, Float64}}
) where R
    # We add a bonus for individuals that sort all test cases correctly
    # The extra logic is to make this work for any record with a fitness field
    new_evaluations = Vector{R}(undef, length(records))
    fields = fieldnames(R)
    fitness_idx = findfirst(field -> field == :fitness, fields)
    for (i, record) in enumerate(records)
        if minimum(values(outcomes[record.id])) == evaluator.num_tests_per_parasite
            vs = [getfield(record, field) for field in fields]
            vs[fitness_idx] = record.fitness + bonuses[i]
            new_evaluations[i] = R(vs...)
        else
            new_evaluations[i] = record
        end
    end
    return new_evaluations
end

function CoEvo.Evaluators.evaluate(
    evaluator::SortingNetworkEvaluator,
    rng::AbstractRNG,
    species::AbstractSpecies,
    outcomes::Dict{Int, SortedDict{Int, Float64}}
) 
    @assert species.population[1].genotype isa SortingNetworkGenotype
    # But we add a bonus for individuals that hit evaluator.max_sort_fitness
    # proportional to how close they are to the minimum number of codons

    # Align the individuals with the records
    evaluation = evaluate(evaluator.evaluator, rng, species, outcomes)
    ind_dict = Dict(ind.id=>ind for ind in [species.population; species.children])
    individuals = [ind_dict[record.id] for record in evaluation.records]

    # Compute the bonuses for each individual
    bonuses = [(ind.genotype.max_codons - length(ind.genotype.codons))/
               (ind.genotype.max_codons - ind.genotype.min_codons)
               for ind in individuals]

    # Add the bonuses to the fitnesses
    new_records = add_bonus_fitness(evaluator, evaluation.records, bonuses, outcomes)

    # Sort the records by new fitnesses
    sort!(new_records, by = x -> x.fitness, rev = true)

    evaluation = OutcomeScalarFitnessEvaluation(species.id, new_records, outcomes)
    return evaluation
end
end
