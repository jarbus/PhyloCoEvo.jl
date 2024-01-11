module EstimatorEcosystem

export EstimatorEcosystemCreator, create_ecosystem, evolve!

import CoEvo.Ecosystems: create_ecosystem, evolve!

using Serialization
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
    check_interval = 1000,
    check_name = "checkpoint.jls"
)
    last_reproduce_time = 0.0

    # TODO This is bad and should be changed
    # This code is just to get the experiment directory
    # for use in checkpointing, we assume the basic archiver
    # follows the path XDIR/data/archiver.jld2, and we save
    # the checkpoint in XDIR/data
    archiver_path = estimator_ecosystem_creator.basic.archiver.archive_path
    data_dir = dirname(archiver_path) 
    @assert endswith(data_dir, "data")

    check_path = joinpath(data_dir, check_name)
    tmp_check_path = check_path*"-tmp"
    archiver_check_path = joinpath(data_dir, "archiver-checkpoint.jld2")
    tmp_archiver_check_path = joinpath(data_dir, "archiver-checkpoint.jld2-tmp")

    if isfile(check_path)
	println("DESERIALIZING CHECKPOINT")
        check = deserialize(check_path)
	estimator_ecosystem_creator = check["ecosystem_creator"]
	ecosystem = check["ecosystem"]
	start_gen = check["generation"] + 1
    	basic = estimator_ecosystem_creator.basic
	mv(archiver_check_path, archiver_path, force=true)
	println("RESUMING FROM $start_gen")
    else
	println("CREATING NEW RUN")
    	basic = estimator_ecosystem_creator.basic
    	ecosystem = create_ecosystem(basic)
	start_gen = 1
    end

    for generation in start_gen:n_generations
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
	if generation % check_interval == 1
	    println("Serializing on generation $generation")
	    # We write to a temporary path, then overwrite the previous checkpoint 
	    # by renaming the file. We are almost guaranteed to never corrupt a prior
	    # checkpoint by renaming the file.
	    serialize(tmp_check_path, Dict(
			          "ecosystem_creator"=>estimator_ecosystem_creator,
				  "ecosystem" => ecosystem,
				  "generation" => generation))
        
	    println("Verifying checkpoint integrity on $generation")
	    Serialization.deserialize(tmp_check_path)     # Check checkpoint integrity 
	    println("Integrity verified")
	    # make a checkpoint of the archiver without overriding previous checkpoint
	    cp(archiver_path, tmp_archiver_check_path, force=true)
	    # quickly rename temporary checkpoints to main checkpoint files
	    mv(tmp_archiver_check_path, archiver_check_path, force=true)
	    mv(tmp_check_path, check_path, force=true)  # If we can deserialize with no errors
				                  # then overwrite prior check
	end
    end

    return ecosystem
end

end
