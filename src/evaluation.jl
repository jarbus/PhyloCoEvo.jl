export OutcomeScalarFitnessEvaluation, OutcomeScalarFitnessEvaluator, OutcomeNSGAIIEvaluation, OutcomeNSGAIIEvaluator

using DataStructures: SortedDict
using CoEvo.Species: AbstractSpecies
using CoEvo.Individuals: Individual
using CoEvo.Evaluators: evaluate
using CoEvo.Evaluators.NSGAII: NSGAIIRecord, NSGAIIEvaluation, NSGAIIEvaluator
using CoEvo.Evaluators.ScalarFitness: ScalarFitnessEvaluator, ScalarFitnessEvaluation, ScalarFitnessRecord
using CoEvo.Observers: Observation

struct OutcomeScalarFitnessEvaluation <: CoEvo.Evaluators.Evaluation
    species_id::String
    records::Vector{ScalarFitnessRecord}
    outcomes::Dict{Int, Dict{Int, Float64}}
end

struct OutcomeNSGAIIEvaluation <: Evaluation
    species_id::String
    records::Vector{NSGAIIRecord}
    outcomes::Dict{Int, Dict{Int, Float64}}
end

Base.@kwdef struct OutcomeScalarFitnessEvaluator <: CoEvo.Evaluators.Evaluator 
    maximize::Bool = true
    epsilon::Float64 = 1e-6
end

Base.@kwdef struct OutcomeNSGAIIEvaluator <: Evaluator 
    maximize::Bool = true
    perform_disco::Bool = true
    max_clusters::Int = -1
    function_minimums::Union{Vector{Float64}, Nothing} = nothing
    function_maximums::Union{Vector{Float64}, Nothing} = nothing
end

function CoEvo.Evaluators.evaluate(
    evaluator::OutcomeScalarFitnessEvaluator,
    rng::AbstractRNG,
    species::AbstractSpecies,
    outcomes::Dict{Int, SortedDict{Int, Float64}}
) 
    """Wrapper around ScalarFitnessEvaluator to create OutcomeScalarFitnessEvaluation
    """
    scalar_fitness_evaluator = ScalarFitnessEvaluator(
        maximize = evaluator.maximize,
        epsilon = evaluator.epsilon
    )
    scalar_fitness_evaluation = evaluate(
        scalar_fitness_evaluator, rng, species, outcomes
    )
    evaluation = OutcomeScalarFitnessEvaluation(species.id, scalar_fitness_evaluation.records, outcomes)
    return evaluation
end

function CoEvo.Evaluators.evaluate(
    evaluator::OutcomeNSGAIIEvaluator,
    random_number_generator::AbstractRNG,
    species::AbstractSpecies,
    outcomes::Dict{Int, SortedDict{Int, Float64}})
    """Wrapper around NSGAIIEvaluator to create OutcomeNSGAIIEvaluation
    """
    nsgaii_evaluator = NSGAIIEvaluator(
        maximize = evaluator.maximize,
        perform_disco = evaluator.perform_disco,
        max_clusters = evaluator.max_clusters,
        function_minimums = evaluator.function_minimums,
        function_maximums = evaluator.function_maximums
    )
    nsgaii_evaluation = evaluate(
        nsgaii_evaluator, random_number_generator, species, outcomes
    )
    evaluation = OutcomeNSGAIIEvaluation(species.id, nsgaii_evaluation.records, outcomes)
    return evaluation
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


function CoEvo.Ecosystems.Basic.evaluate_species(
    evaluators::Vector{<:Evaluator},
    random_number_generator::AbstractRNG,
    species::Vector{<:AbstractSpecies},
    individual_outcomes::Dict{Int, SortedDict{Int, Float64}},
    observations::Vector{<:Observation},
)
    estimate_outcomes!(individual_outcomes, evaluators, species, observations)

    evaluations = [
        evaluate(evaluator, random_number_generator, species, individual_outcomes)
        for (evaluator, species) in zip(evaluators, species)
    ]

    return evaluations
end



function CoEvo.Replacers.replace(
    replacer::CoEvo.Replacers.Generational.GenerationalReplacer,
    rng::AbstractRNG, 
    species::AbstractSpecies,
    evaluation::OutcomeScalarFitnessEvaluation
)
    println("len evaluation.records: ", length(evaluation.records))
    sf_evaluation = ScalarFitnessEvaluation(evaluation.species_id, evaluation.records)
    new_population = CoEvo.Replacers.replace(replacer, rng, species, sf_evaluation)
    return new_population
end

function CoEvo.Selectors.select(
    selector::CoEvo.Selectors.FitnessProportionate.FitnessProportionateSelector,
    random_number_generator::AbstractRNG, 
    new_population::Vector{<:Individual},
    evaluation::OutcomeScalarFitnessEvaluation
)
    sf_evaluation = ScalarFitnessEvaluation(evaluation.species_id, evaluation.records)
    selected_population = CoEvo.Selectors.select(selector, random_number_generator, new_population, sf_evaluation)
    return selected_population
end
