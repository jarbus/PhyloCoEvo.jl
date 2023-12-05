module Vectors

export TwoIndexUniformVectorMutator

import CoEvo.Mutators: mutate

using Random: AbstractRNG, randn
using CoEvo.Counters: Counter
using CoEvo.Genotypes.Vectors: BasicVectorGenotype
using CoEvo.Mutators: Mutator
using ...Utils: two_rand

Base.@kwdef struct TwoIndexUniformVectorMutator <: Mutator
    lower_bound::Float64 = -0.1
    upper_bound::Float64 = 0.1
end

function mutate(
    mutator::TwoIndexUniformVectorMutator, 
    random_number_generator::AbstractRNG,
    ::Counter,
    genotype::BasicVectorGenotype{T}
) where T

    # select two random indices
    idxs = two_rand(random_number_generator, 1:length(genotype)) |> collect
    # generate noise vector
    noise_vector = rand(random_number_generator, length(idxs))
    # rescale noise vector to be in range [lower_bound, upper_bound]
    noise_vector = noise_vector .* (mutator.upper_bound - mutator.lower_bound) .+ mutator.lower_bound
    # mutate genotype
    mutated_genes = copy(genotype.genes)
    mutated_genes[idxs] = mutated_genes[idxs] .+ noise_vector
    mutated_genotype = BasicVectorGenotype(mutated_genes)
    return mutated_genotype
end

end

