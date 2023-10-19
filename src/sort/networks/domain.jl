# TODO fill this in, can use the stateless environment creator
using CoEvo.Ecosystems.Metrics.Outcomes.Types.NumbersGame: NumbersGame as NumbersGameMetrics
using .NumbersGameMetrics: NumbersGameMetric, Control, Sum, Gradient, Focusing, Relativism
using CoEvo.Ecosystems.Interactions.Domains.Abstract: Domain
using CoEvo.Ecosystems.Metrics.Outcomes.Abstract: OutcomeMetric

abstract type AbstractSortingNetworkMetric <: OutcomeMetric end
# We can add more metrics here
struct SortingNetworkMetric <: AbstractSortingNetworkMetric end

struct SortingNetworkDomain{O <: AbstractSortingNetworkMetric} <: Domain{O}
    outcome_metric::O
end
