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
        outcome_set = CoEvo.Environments.get_outcome_set(environment.domain.outcome_metric,
                                                         result,
                                                         environment.phenotypes[1])
        return outcome_set
    elseif environment.phenotypes[2] isa SortingNetworkPhenotype
        result = netsort(environment.phenotypes[2], environment.phenotypes[1])
        outcome_set = CoEvo.Environments.get_outcome_set(environment.domain.outcome_metric,
                                                         result,
                                                         environment.phenotypes[1])
        return outcome_set[2:-1:1]
    else
        error("Neither phenotype is a sorting network phenotype")
    end
end

function CoEvo.Environments.get_outcome_set(::Partial, results::Vector{Vector{Int64}}, snp::SortingNetworkPhenotype)
    # For partial, we just want the number of results sorted correctly
    num_correct = 0
    for r in results
        correct = true
        for i in 1:length(r)-1
            if r[i] > r[i+1]
                correct = false
                break
            end
        end
        if correct
            num_correct += 1
        end
    end
    # If all results are sorted correctly, we add a bonus for how close we are
    # to the minimum number of codons
    if num_correct == length(results)
        if snp.max_codons == snp.min_codons
            percent_extra_swaps = 0
        else
            percent_extra_swaps = ((size(snp.network, 1) - snp.min_codons) / (snp.max_codons - snp.min_codons))
        end
        return Float64[num_correct + (1-percent_extra_swaps),
                       length(results) - num_correct + percent_extra_swaps]
    else
        return Float64[num_correct, length(results) - num_correct]
    end

end
