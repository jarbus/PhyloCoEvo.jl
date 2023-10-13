export OutcomeScalarFitnessEvaluation, OutcomeScalarFitnessEvaluator

using DataStructures: OrderedDict
using CoEvo.Species.Abstract: AbstractSpecies
using CoEvo.Species.Individuals: Individual
using CoEvo.Species.Evaluators.Abstract: Evaluation, Evaluator
using CoEvo.Species.Evaluators.Interfaces: create_evaluation, get_ranked_ids

struct OutcomeScalarFitnessEvaluation <: Evaluation
    species_id::String
    fitnesses::OrderedDict{Int, Float64}
    outcomes::Dict{Int, Dict{Int, Float64}}
end

Base.@kwdef struct OutcomeScalarFitnessEvaluator <: Evaluator 
    maximize::Bool = true
    epsilon::Float64 = 1e-6
end

function CoEvo.create_evaluation(
    evaluator::OutcomeScalarFitnessEvaluator,
    species::AbstractSpecies,
    outcomes::Dict{Int, Dict{Int, Float64}}
) 
    """Wrapper around ScalarFitnessEvaluator to create OutcomeScalarFitnessEvaluation
    """
    scalar_fitness_evaluator = ScalarFitnessEvaluator(
        maximize = evaluator.maximize,
        epsilon = evaluator.epsilon
    )
    scalar_fitness_evaluation = create_evaluation(
        scalar_fitness_evaluator, species, outcomes
    )
    evaluation = OutcomeScalarFitnessEvaluation(species.id, scalar_fitness_evaluation.fitnesses, outcomes)
    return evaluation
end

function CoEvo.select(
    selector::FitnessProportionateSelector,
    rng::AbstractRNG, 
    new_pop::Dict{Int, <:Individual},
    evaluation::OutcomeScalarFitnessEvaluation
)
    """Wrapper around FitnessProportionateSelector to select individuals given OutcomeScalarFitnessEvaluation
    """
    scalar_fitness_evaluation = ScalarFitnessEvaluation(
        evaluation.species_id,
        evaluation.fitnesses,
        []
    )
    return select(selector, rng, new_pop, scalar_fitness_evaluation)
end

function CoEvo.replace(
    replacer::GenerationalReplacer,
    rng::AbstractRNG, 
    species::AbstractSpecies,
    evaluation::OutcomeScalarFitnessEvaluation
)
    """Wrapper around GenerationalReplacer to replace individuals given OutcomeScalarFitnessEvaluation
    """
    scalar_fitness_evaluation = ScalarFitnessEvaluation(evaluation.species_id, evaluation.fitnesses, [])
    return replace(replacer, rng, species, scalar_fitness_evaluation)
end


function CoEvo.get_ranked_ids(evaluator::OutcomeScalarFitnessEvaluation, ids::Vector{Int})
    ranked_ids = filter(
        indiv_id -> indiv_id in ids, keys(evaluator.fitnesses)
    )
    return ranked_ids
end

