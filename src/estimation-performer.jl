using Distributed: remotecall, fetch
using CoEvo.Interactions: Interaction, interact
using CoEvo.Performers: Performer, perform
using CoEvo.Jobs: Job
using CoEvo.Jobs.Basic: BasicJob
using CoEvo.Phenotypes: Phenotype
using CoEvo.Environments: create_environment 
using CoEvo.Results: Result
using CoEvo.Observers: Observer

export EstimationPerformer

Base.@kwdef struct EstimationPerformer <: Performer 
    n_workers::Int
end

function CoEvo.Performers.perform(::EstimationPerformer, job::BasicJob)
    results = Result[]
    for match in job.matches
        interaction = job.interactions[match.interaction_id]
        phenotypes = Phenotype[
            job.phenotypes[individual_id] for individual_id in match.individual_ids
        ]
        result = interact(
            interaction,
            match.individual_ids,
            phenotypes
        )
        push!(results, result)
    end
    return results
end

function CoEvo.Performers.perform(performer::EstimationPerformer, jobs::Vector{<:BasicJob})
    if length(jobs) == 1
        results = perform(performer, jobs[1])
    else
        futures = [remotecall(perform, i, performer, job) for (i, job) in enumerate(jobs)]
        results = [fetch(f) for f in futures]
    end
    results::Vector{Result} = vcat(results...)
    return results
end
