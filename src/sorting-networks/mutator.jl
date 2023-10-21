export SortingNetworkMutator, SortingNetworkTestCaseMutator, mutate
using Random: AbstractRNG, randn

using CoEvo


Base.@kwdef struct SortingNetworkMutator <: Mutator 
    activation_prob::Float64 = 0.1
    bit_flip_prob::Float64 = 0.01
end

Base.@kwdef struct SortingNetworkTestCaseMutator <: Mutator 
    swap_prob::Float64 = 0.1
end

function CoEvo.mutate(
    mutator::SortingNetworkMutator, 
    rng::AbstractRNG, 
    ::Counter,
    geno::SortingNetworkGenotype
)
    # TODO implement SortingNetworkMutator
end

function CoEvo.mutate(
    mutator::SortingNetworkTestCaseMutator, 
    rng::AbstractRNG, 
    ::Counter,
    geno::SortingNetworkTestCaseGenotype
)
    # TODO implement SortingNetworkTestCaseMutator
end
