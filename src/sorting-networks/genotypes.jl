export SortingNetworkGenotypeCreator, SortingNetworkGenotype, create_genotypes

using Random: AbstractRNG
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
