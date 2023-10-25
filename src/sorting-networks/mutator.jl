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

function CoEvo.Mutators.mutate(
    ::SortingNetworkMutator, 
    rng::AbstractRNG, 
    gene_id_counter::Counter,
    geno::SortingNetworkGenotype
)
    # Perform one activation flip
    new_codons = [codon for codon in geno.codons]
    index = rand(rng, 1:length(new_codons))
    activation_id = new_codons[index].id
    activation_data = new_codons[index].data ⊻ (1 << rand(rng, 7:15))
    new_codons[index] = SortingNetworkCodon(activation_id, activation_data)
    # Perform one bit flip
    active_codons = is_active.(new_codons)
    if any(active_codons)
        index = rand(rng, findall(active_codons))
        bitflip_data = new_codons[index].data ⊻ (1 << rand(rng, 0:7))
        bitflip_id = count!(gene_id_counter)
        new_codons[index] = SortingNetworkCodon(bitflip_id, bitflip_data)
    end
    return SortingNetworkGenotype(Tuple(new_codons), geno.n_inputs)
end

function CoEvo.Mutators.mutate(
    ::SortingNetworkTestCaseMutator, 
    rng::AbstractRNG, 
    ::Counter,
    geno::SortingNetworkTestCaseGenotype
)
    # swap two inputs
    new_inputs = [i for i in geno.inputs]
    index = rand(rng, 1:length(new_inputs))
    swap_index = rand(rng, setdiff(1:length(new_inputs), index))
    new_inputs[index], new_inputs[swap_index] = new_inputs[swap_index], new_inputs[index]
    return SortingNetworkTestCaseGenotype(geno.id, Tuple(new_inputs))
end
