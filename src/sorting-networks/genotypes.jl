export SortingNetworkGenotypeCreator, SortingNetworkGenotype, create_genotypes, SortingNetworkTestCaseGenotypeCreator, SortingNetworkTestCaseGenotype

using Random: AbstractRNG, shuffle
using CoEvo.Ecosystems.Species.Genotypes.Abstract: Gene
using CoEvo.Ecosystems.Utilities.Counters: next!

struct SortingNetworkCodon <: Gene
    id::Int
    data::UInt16
end

struct SortingNetworkGenotype{N} <: Genotype
    codons::NTuple{N, SortingNetworkCodon}
    n_inputs::Int # number of inputs
end

struct SortingNetworkGenotypeCreator <: GenotypeCreator
    n_codons::Int
    n_inputs::Int
end

struct SortingNetworkTestCaseGenotype{N} <: Genotype
    id::Int
    inputs::NTuple{N, Int64}
end

struct SortingNetworkTestCaseGenotypeCreator <: GenotypeCreator
    n_inputs::Int
end

function is_active(codon::SortingNetworkCodon)
    count_ones(codon.data & 0b1111111100000000) > 3
end
function num_active(geno::SortingNetworkGenotype)
    sum(is_active.(geno.codons))
end

function CoEvo.create_genotypes(
    genotype_creator::SortingNetworkGenotypeCreator,
    rng::AbstractRNG,
    gene_id_counter::Counter,
    n_pop::Int
)
    # TODO: initialize networks with the first half of the codons of successful networks per hillis
    genotypes = [
        SortingNetworkGenotype(
            [ SortingNetworkCodon(next!(gene_id_counter), rand(rng, UInt16))
                for i in 1:genotype_creator.n_codons] |> Tuple,
            genotype_creator.n_inputs
        ) for i in 1:n_pop
    ]

    return genotypes
end

function swap(rng::AbstractRNG, num_inputs::Int, num_swaps::Int)
    inputs = collect(1:num_inputs)
    for i in 1:num_swaps
        # choose two random input indicies
        a = rand(rng, 1:num_inputs)
        b = rand(rng, setdiff(1:num_inputs, a))
        inputs[a], inputs[b] = inputs[b], inputs[a]
    end
    return inputs
end

function CoEvo.create_genotypes(
        genotype_creator::SortingNetworkTestCaseGenotypeCreator,
        rng::AbstractRNG,
        gene_id_counter::Counter,
        n_pop::Int
    )
    genotypes = [
        SortingNetworkTestCaseGenotype(
            next!(gene_id_counter),
            swap(rng, genotype_creator.n_inputs, 2) |> Tuple
        ) for _ in 1:n_pop ]

    return genotypes
end
