export SortingNetworkGenotypeCreator, SortingNetworkGenotype, create_genotypes

using Random: AbstractRNG
using CoEvo.Ecosystems.Species.Genotypes.Abstract: Gene
using CoEvo.Ecosystems.Utilities.Counters: next!

struct SortingNetworkCodon <: Gene
    id::Int
    bits::UInt8
end

struct SortingNetworkChromosome{N} <: Gene
    codons::NTuple{N, SortingNetworkCodon}
end

struct SortingNetworkGenotype{N} <: Genotype
    chromosomes::NTuple{N, SortingNetworkChromosome}
    n_inputs::Int # number of inputs
end


struct SortingNetworkGenotypeCreator <: GenotypeCreator
    n_chromosomes::Int
    n_codons_per_chromosome::Int
    n_inputs::Int
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
            [ SortingNetworkChromosome(
                [ SortingNetworkCodon(next!(gene_id_counter), rand(rng, UInt8))
                for i in 1:genotype_creator.n_codons_per_chromosome] |> Tuple
            ) for i in 1:genotype_creator.n_chromosomes] |> Tuple,
            genotype_creator.n_inputs
        ) for i in 1:n_pop
    ]

    return genotypes
end
