export SortingNetworkGenotypeCreator, SortingNetworkGenotype, create_genotypes, SortingNetworkTestCaseGenotypeCreator, SortingNetworkTestCaseGenotype

using Random: AbstractRNG, shuffle
using CoEvo.Genotypes: Gene
using CoEvo.Counters: count!
using CoEvo.Genotypes: Genotype

struct SortingNetworkCodon <: Gene
    id::Int
    one::Int
    two::Int
end

struct SortingNetworkGenotype <: Genotype
    codons::Vector{SortingNetworkCodon}
    n_inputs::Int # number of inputs
end

struct SortingNetworkGenotypeCreator <: CoEvo.Genotypes.GenotypeCreator
    n_codons::Int
    n_inputs::Int
end

struct SortingNetworkTestCaseGenotype <: Genotype
    id::Int
    tests::Vector{Vector{Int64}}
end

struct SortingNetworkTestCaseGenotypeCreator <: CoEvo.Genotypes.GenotypeCreator
    n_tests::Int
    n_inputs::Int
end

function two_random_inputs(rng::AbstractRNG, n_inputs::Int)
    input_1 = rand(rng, 1:n_inputs)
    input_2 = input_1
    while true
        input_2 = rand(rng, 1:n_inputs)
        if input_1 != input_2
            break
        end
    end
    return input_1, input_2
end

function random_codon(rng::AbstractRNG,
                      gene_id_counter::Counter,
                      n_inputs::Int)
    id = count!(gene_id_counter)
    input_1, input_2 = two_random_inputs(rng, n_inputs)
    return SortingNetworkCodon(id, input_1, input_2)
end

function CoEvo.Genotypes.create_genotypes(
    genotype_creator::SortingNetworkGenotypeCreator,
    rng::AbstractRNG,
    gene_id_counter::Counter,
    n_pop::Int
)
    # TODO: initialize networks with the first half of the codons of successful networks per hillis
    
    genotypes = [
        SortingNetworkGenotype(
            [random_codon(rng, gene_id_counter, genotype_creator.n_inputs)
                for i in 1:genotype_creator.n_codons],
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

function CoEvo.Genotypes.create_genotypes(
        genotype_creator::SortingNetworkTestCaseGenotypeCreator,
        rng::AbstractRNG,
        gene_id_counter::Counter,
        n_pop::Int
    )
    genotypes = [
        SortingNetworkTestCaseGenotype(
            count!(gene_id_counter),
            [ swap(rng, genotype_creator.n_inputs, 2)
             for i in 1:genotype_creator.n_tests ]
        ) for _ in 1:n_pop ]

    return genotypes
end
