export OutcomeScalarFitnessEvaluation, OutcomeScalarFitnessEvaluator

using DataStructures: OrderedDict
using CoEvo.Species: AbstractSpecies
using CoEvo.Individuals: Individual
using CoEvo.Evaluators: create_evaluation
using CoEvo.Evaluators.ScalarFitness: ScalarFitnessEvaluator, ScalarFitnessEvaluation

struct OutcomeScalarFitnessEvaluation <: CoEvo.Evaluators.Evaluation
    species_id::String
    fitnesses::OrderedDict{Int, Float64}
    outcomes::Dict{Int, Dict{Int, Float64}}
end

Base.@kwdef struct OutcomeScalarFitnessEvaluator <: CoEvo.Evaluators.Evaluator 
    maximize::Bool = true
    epsilon::Float64 = 1e-6
end

function CoEvo.Evaluators.create_evaluation(
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

function CoEvo.Selectors.select(
    selector::CoEvo.Selectors.FitnessProportionate.FitnessProportionateSelector,
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

function CoEvo.Replacers.replace(
    replacer::CoEvo.Replacers.GenerationalReplacer,
    rng::AbstractRNG, 
    species::AbstractSpecies,
    evaluation::OutcomeScalarFitnessEvaluation
)
    """Wrapper around GenerationalReplacer to replace individuals given OutcomeScalarFitnessEvaluation
    """
    scalar_fitness_evaluation = ScalarFitnessEvaluation(evaluation.species_id, evaluation.fitnesses, [])
    return replace(replacer, rng, species, scalar_fitness_evaluation)
end
