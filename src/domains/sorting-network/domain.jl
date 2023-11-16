module SortingNetwork
export SortingNetworkDomain, Partial, get_outcome_set, PartialPlusBonus
using CoEvo
using CoEvo.Domains: Domain
using CoEvo.Metrics: Metric
using CoEvo.Phenotypes: Phenotype
using ...Phenotypes.SortingNetwork: SortingNetworkPhenotype, netsort
import CoEvo.Environments.get_outcome_set

abstract type AbstractSortingNetworkMetric <: Metric end
# We can add more metrics here
struct Partial <: AbstractSortingNetworkMetric
"""Give the number of results sorted correctly and incorrectly"""
end

struct PartialPlusBonus <: AbstractSortingNetworkMetric 
"""Give the number of results sorted correctly and incorrectly. If sorted correctly, add extra reward based on how small the network is."""
end


struct SortingNetworkDomain{O <: AbstractSortingNetworkMetric} <: CoEvo.Domains.Domain{O}
    outcome_metric::O
end

function get_outcome_set(
        environment::CoEvo.Environments.Stateless.StatelessEnvironment{D, Phenotype}) where {D <: SortingNetworkDomain}

    if environment.phenotypes[1] isa SortingNetworkPhenotype
        result = netsort(environment.phenotypes[1], environment.phenotypes[2])
        outcome_set = get_outcome_set(environment.domain.outcome_metric,
                                                         result,
                                                         environment.phenotypes[1])
        return outcome_set
    elseif environment.phenotypes[2] isa SortingNetworkPhenotype
        result = netsort(environment.phenotypes[2], environment.phenotypes[1])
        outcome_set = get_outcome_set(environment.domain.outcome_metric,
                                                         result,
                                                         environment.phenotypes[2])
        return outcome_set[2:-1:1]
    else
        error("Neither phenotype is a sorting network phenotype")
    end
end

function count_correct(results::Vector{Vector{Int64}})
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
    num_correct
end

function get_outcome_set(::Partial, results::Vector{Vector{Int64}}, snp::SortingNetworkPhenotype)
    # For partial, we just want the number of results sorted correctly
    num_correct = count_correct(results)
    Float64[num_correct, length(results) - num_correct]
end

function get_outcome_set(::PartialPlusBonus, results::Vector{Vector{Int64}}, snp::SortingNetworkPhenotype)
    # For partial, we just want the number of results sorted correctly
    num_correct = count_correct(results) 
    if num_correct != length(results)
        return Float64[num_correct, length(results) - num_correct]
    end
    sn_size = size(snp.network, 1)
    bonus = (snp.max_codons - sn_size) / (snp.max_codons - snp.min_codons)
    @assert 0 <= bonus <= 1 "Bonus should be between 0 and 1, but is $bonus"
    return Float64[num_correct + bonus, 1 - bonus]
end

end
