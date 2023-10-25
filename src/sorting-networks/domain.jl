export SortingNetworkDomain, Partial
using CoEvo.Domains: Domain
using CoEvo.Metrics: Metric

abstract type AbstractSortingNetworkMetric <: Metric end
# We can add more metrics here
struct Partial <: AbstractSortingNetworkMetric end

struct SortingNetworkDomain{O <: AbstractSortingNetworkMetric} <: CoEvo.Domains.Domain{O}
    outcome_metric::O
end

function CoEvo.Environments.get_outcome_set(
        environment::CoEvo.Environments.Stateless.StatelessEnvironment{D, Phenotype}) where {D <: SortingNetworkDomain}

    if environment.phenotypes[1] isa SortingNetworkPhenotype
        result = netsort(environment.phenotypes[1], environment.phenotypes[2])
        outcome_set = get_outcome_set(environment.domain.outcome_metric, result)
        return outcome_set
    elseif environment.phenotypes[2] isa SortingNetworkPhenotype
        result = netsort(environment.phenotypes[2], environment.phenotypes[1])
        outcome_set = get_outcome_set(environment.domain.outcome_metric, result)
        return outcome_set[2:-1:1]
    else
        error("Neither phenotype is a sorting network phenotype")
    end
end

function CoEvo.Environments.get_outcome_set(::Partial, results::Vector{<:Int64})
    # For partial, we just want the number of correct numbers
    correct = sum(results .== 1:length(results))
    return Float64[correct, length(results) - correct]
end
