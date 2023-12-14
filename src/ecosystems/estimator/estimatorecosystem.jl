module EstimatorEcosystem

export EstimatorEcosystemCreator, create_ecosystem, evolve!

import CoEvo.Ecosystems: create_ecosystem, evolve!

using DataStructures: SortedDict
using Random: AbstractRNG
using StableRNGs: StableRNG
using CoEvo.Species: AbstractSpecies
using CoEvo.Counters: Counter, Basic.BasicCounter
using CoEvo.Evaluators: Evaluation, Evaluator, evaluate
using CoEvo.SpeciesCreators: SpeciesCreator, create_species
using CoEvo.SpeciesCreators.Basic: BasicSpeciesCreator
using CoEvo.Jobs: JobCreator, create_jobs
using CoEvo.Performers: Performer
using CoEvo.Interactions: Interaction
using CoEvo.Results: Result, get_individual_outcomes, get_observations
using CoEvo.Observers: Observation
using CoEvo.Observers.Null: NullObservation
using CoEvo.Reporters.Runtime: RuntimeReporter, create_runtime_report
using CoEvo.Reporters: Reporter, Report, create_report
using CoEvo.Archivers: Archiver, archive!, archive_reports!
using CoEvo.Performers: perform
using CoEvo.States.Basic: BasicCoevolutionaryStateCreator, BasicCoevolutionaryState
using CoEvo.States: State, StateCreator
using CoEvo.Ecosystems: Ecosystem, EcosystemCreator
using CoEvo.Ecosystems.Basic: BasicEcosystemCreator, BasicEcosystem
using ...Estimators: Estimator, estimate!
using CoEvo.Ecosystems.Basic: evaluate_species, construct_new_species, create_state, create_all_reports
import CoEvo.Ecosystems: evolve!, create_ecosystem


Base.@kwdef struct EstimatorEcosystemCreator <: EcosystemCreator
    basic::EcosystemCreator
    estimators::Vector{<:Estimator}
end

show(io::IO, c::EstimatorEcosystemCreator) = show(io, c.basic)
# Forwards all properties to basic except for estimators
Base.getproperty(e::EstimatorEcosystemCreator, p::Symbol) =
    p ∈ (:estimators, :basic) ?  getfield(e,p) : getproperty(getfield(e,:basic),p)

   
function create_ecosystem(
    ecosystem_creator::EstimatorEcosystemCreator,
    gen::Int, 
    ecosystem::Ecosystem, 
    results::Vector{<:Result}, 
    reports::Vector{Report}
)
    individual_outcomes = get_individual_outcomes(results)

    start = time()
    estimate!(ecosystem_creator.estimators, individual_outcomes, ecosystem.species)
    println("Estimation time: ", time() - start)

    observations = Observation[]

    evaluations = evaluate_species(
        ecosystem_creator.basic, ecosystem, individual_outcomes, observations
    )

    state = create_state(
        ecosystem_creator.basic.state_creator, 
        ecosystem_creator.basic, 
        gen, 
        ecosystem, 
        individual_outcomes,
        evaluations,
        observations,
    )
    generation_reports = create_all_reports(state, ecosystem_creator.basic.reporters)
    append!(reports, generation_reports)
    archive_reports!(ecosystem_creator.basic.archiver, reports)
    if gen % ecosystem_creator.basic.garbage_collection_interval == 0
        Base.GC.gc()
    end
    all_new_species = construct_new_species(state, ecosystem_creator.basic.species_creators)
    new_eco = BasicEcosystem(ecosystem_creator.basic.id, all_new_species)
    
    return new_eco
end
using JLD2: @save


function evolve!(
    estimator_ecosystem_creator::EstimatorEcosystemCreator;
    n_generations::Int = 100,
)
    basic = estimator_ecosystem_creator.basic
    ecosystem = create_ecosystem(estimator_ecosystem_creator.basic)
    last_reproduce_time = 0.0
    for generation in 1:n_generations
        eval_time_start = time()
        phenotype_creators = [
            species_creator.phenotype_creator 
            for species_creator in basic.species_creators
        ]
        jobs = create_jobs(
            basic.job_creator,
            basic.random_number_generator, 
            ecosystem.species,
            phenotype_creators,
        )
        results = perform(basic.performer, jobs)
        eval_time = time() - eval_time_start
        runtime_report = create_runtime_report(
            basic.runtime_reporter, 
            basic.id, 
            generation, 
            eval_time, 
            last_reproduce_time
        )
        reports = Report[runtime_report]

        last_reproduce_time_start = time()
        ecosystem = create_ecosystem(estimator_ecosystem_creator, generation, ecosystem, results, reports)
        last_reproduce_time = time() - last_reproduce_time_start
    end

    return ecosystem
end

end
