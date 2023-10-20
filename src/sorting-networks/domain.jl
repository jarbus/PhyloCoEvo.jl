export SortingNetworkDomain, Partial
using CoEvo.Ecosystems.Metrics.Outcomes.Types.NumbersGame: NumbersGame as NumbersGameMetrics
using .NumbersGameMetrics: NumbersGameMetric, Control, Sum, Gradient, Focusing, Relativism
using CoEvo.Ecosystems.Interactions.Domains.Abstract: Domain
using CoEvo.Ecosystems.Metrics.Outcomes.Abstract: OutcomeMetric

abstract type AbstractSortingNetworkMetric <: OutcomeMetric end
# We can add more metrics here
struct Partial <: AbstractSortingNetworkMetric end

struct SortingNetworkDomain{O <: AbstractSortingNetworkMetric} <: Domain{O}
    outcome_metric::O
end

function get_outcome_set(
    environment::StatelessEnvironment{D, <:SortingNetworkPhenotype}) where {D <: SortingNetworkDomain}

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

function get_outcome_set(::Partial, results::Vector{<:Int64})
    # For partial, we just want the number of correct numbers
    return sum(results .== 1:length(results))
end
