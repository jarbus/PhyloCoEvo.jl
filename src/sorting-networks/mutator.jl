export SortingNetworkMutator, SortingNetworkTestCaseMutator, mutate
using Random: AbstractRNG, randn

using CoEvo


Base.@kwdef struct SortingNetworkMutator <: Mutator 
    num_ops_per_mut::Int = 1
end

Base.@kwdef struct SortingNetworkTestCaseMutator <: Mutator 
    num_swaps_per_mut::Int = 1
end

function CoEvo.Mutators.mutate(
    mutator::SortingNetworkMutator, 
    rng::AbstractRNG, 
    gene_id_counter::Counter,
    geno::SortingNetworkGenotype
)
    new_codons = [c for c in geno.codons]
    for _ in 1:mutator.num_ops_per_mut
        choice = rand(rng, 1:4)
        if choice == 1
            # insert
            length(new_codons) >= geno.max_codons && continue
            new_codon = random_codon(rng, gene_id_counter, geno.n_inputs)
            random_insert!(new_codons, new_codon)
        elseif choice == 2
            # delete
            length(new_codons) <= geno.min_codons && continue
            index = rand(rng, 1:length(new_codons))
            deleteat!(new_codons, index)
        elseif choice == 3
            # move
            index = rand(rng, 1:length(new_codons))
            new_index = rand(rng, 1:length(new_codons))
            codon = popat!(new_codons, index)
            insert!(new_codons, new_index, codon)
        elseif choice == 4
            # rewire 
            index = rand(rng, 1:length(new_codons))
            new_id = count!(gene_id_counter)
            new_one, new_two = two_rand(rng, 1:geno.n_inputs)
            new_codons[index] = SortingNetworkCodon(new_id, new_one, new_two)
        end
    end
    SortingNetworkGenotype(new_codons, geno.n_inputs, geno.min_codons, geno.max_codons)
end

function random_insert!(codons::Vector{SortingNetworkCodon}, codon::SortingNetworkCodon)
    index = rand(1:length(codons))
    insert!(codons, index, codon)
end

function CoEvo.Mutators.mutate(
    mut::SortingNetworkTestCaseMutator, 
    rng::AbstractRNG, 
    ::Counter,
    geno::SortingNetworkTestCaseGenotype
)
    # swap two inputs
    new_tests = [[i for i in t] for t in geno.tests]
    for _ in 1:mut.num_swaps_per_mut 
        test_index = rand(rng, 1:length(new_tests))
        swap_index1 = rand(rng, 1:length(new_tests[test_index]))
        swap_index2 = rand(rng, setdiff(1:length(new_tests[test_index]), swap_index1))
        new_tests[test_index][swap_index2], new_tests[test_index][swap_index1] = 
            new_tests[test_index][swap_index1], new_tests[test_index][swap_index2]
    end
    return SortingNetworkTestCaseGenotype(geno.id, new_tests)
end
