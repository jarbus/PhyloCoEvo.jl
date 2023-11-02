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

function create_random_sorting_network(rng::AbstractRNG,
                                       gene_id_counter::Counter,
                                       n_codons::Int,
                                       n_inputs::Int)
    codons = [
        random_codon(rng, gene_id_counter, n_inputs)
        for _ in 1:n_codons
    ]
    return SortingNetworkGenotype(codons, n_inputs)
end

function create_seeded_sorting_network_16(rng::AbstractRNG,
                                          gene_id_counter::Counter,
                                          n_codons::Int)
    # Seed the network with the first half of the codons of the successful
    # networks from Hillis' paper
    swaps = [ 
        # Column 1: swap adjacent
        1  2; 3  4; 5  6; 7  8; 9 10; 11 12; 13 14; 15 16;
        # Columns 2, 3: swap 1 apart
        1 3; 5 7; 9 11; 13 15;
        2 4; 6 8; 10 12; 14 16;
        # Columns 4-7: Swap 3 apart
        1 5; 9 13;
        2 6; 10 14;
        3 7; 11 15;
        4 8; 12 16;
        # Columns 8-15: Swap 7 apart
        1 9; 2 10; 3 11; 4 12; 5 13; 6 14; 7 15; 8 16;
    ]
    codons = [
        SortingNetworkCodon(count!(gene_id_counter), swap[1], swap[2])
        for swap in eachrow(swaps)
    ]
    if n_codons > length(codons)
        codons = vcat(codons, [
            random_codon(rng, gene_id_counter, 16)
            for i in 1:(n_codons - length(codons))
        ])
    else
        error("Not enough codons to seed the network: ncodons: $n_codons < num_seed_codons: $(length(codons))")
    end
    return SortingNetworkGenotype(codons, 16)
end


function CoEvo.Genotypes.create_genotypes(
    genotype_creator::SortingNetworkGenotypeCreator,
    rng::AbstractRNG,
    gene_id_counter::Counter,
    n_pop::Int
)
    # Seed the network with the first half of the codons of the successful
    # networks from Hillis' paper if there are 16 inputs
    if genotype_creator.n_inputs == 16
        genotypes = [ create_seeded_sorting_network_16(
                rng, gene_id_counter, genotype_creator.n_codons)
             for i in 1:n_pop]
        return genotypes
    else
        genotypes = [ create_random_sorting_network(
            rng, gene_id_counter, genotype_creator.n_codons, genotype_creator.n_inputs)
        for i in 1:n_pop]
        return genotypes
    end

end

function swap(rng::AbstractRNG, num_inputs::Int, num_swaps::Int)
    inputs = collect(1:num_inputs)
    for i in 1:num_swaps
        # choose two random input indicies
        a, b = two_random_inputs(rng, num_inputs)
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
