using Distributed: remotecall, fetch
using CoEvo.Ecosystems.Interactions.Results: Result

using CoEvo.Ecosystems.Performers.Abstract: Performer
using CoEvo.Ecosystems.Jobs.Abstract: Job
using CoEvo.Ecosystems.Jobs.Basic: BasicJob
using CoEvo.Ecosystems.Interactions.Abstract: Interaction
using CoEvo.Ecosystems.Species.Phenotypes.Abstract: Phenotype
using CoEvo.Ecosystems.Interactions.Observers.Abstract: Observer
using CoEvo.Ecosystems.Interactions.Environments.Interfaces: create_environment, create_observer
using CoEvo.Ecosystems.Interactions.Methods.Interact: interact

import CoEvo.Ecosystems.Performers.Interfaces: perform

Base.@kwdef struct EstimationPerformer <: Performer 
    n_workers::Int
end

"""
    perform(job::BasicJob) -> Vector{InteractionResult}

Execute the given `job`, which contains various interaction recipes. Each recipe denotes 
specific entities to interact in a domain. The function processes these interactions and 
returns a list of their results.

# Arguments
- `job::BasicJob`: The job containing details about the interactions to be performed.

# Returns
- A `Vector` of `InteractionResult` instances, each detailing the outcome of an interaction.
"""
function CoEvo.perform(::EstimationPerformer, job::BasicJob)
    results = Result[]
    for match in job.matches
        environment_creator = job.interactions[match.interaction_id].environment_creator
        observer_creators = job.interactions[match.interaction_id].observers
        phenotypes = [job.phenotypes[indiv_id] for indiv_id in match.indiv_ids]
        environment = create_environment(environment_creator, phenotypes)
        observers = [create_observer(creator, environment) for creator in observer_creators]
        result = interact(
            match.interaction_id, 
            match.indiv_ids,
            environment, 
            observers
        )
        push!(results, result)
    end
    return results
end


function CoEvo.perform(performer::EstimationPerformer, jobs::Vector{<:Job})
    if length(jobs) == 1
        results = perform(performer, jobs[1])
    else
        futures = [remotecall(perform, i, performer, job) for (i, job) in enumerate(jobs)]
        results = [fetch(f) for f in futures]
    end
    results = vcat(results...)
    return results
end
