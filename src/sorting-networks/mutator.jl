export SortingNetworkMutator, SortingNetworkTestCaseMutator, mutate
using Random: AbstractRNG, randn

using CoEvo


Base.@kwdef struct SortingNetworkMutator <: Mutator 
    min_codons::Int = 60
    max_codons::Int = 80
    num_swaps_per_mut::Int = 1
    num_insert_delete_move_per_mut::Int = 1
end

Base.@kwdef struct SortingNetworkTestCaseMutator <: Mutator 
    swap_prob::Float64 = 0.1
end

function CoEvo.Mutators.mutate(
    mutator::SortingNetworkMutator, 
    rng::AbstractRNG, 
    gene_id_counter::Counter,
    geno::SortingNetworkGenotype
)
    new_codons = [c for c in geno.codons]
    # insert/delete
    for _ in 1:mutator.num_insert_delete_move_per_mut
        choice = rand(rng, 1:3)
        if choice == 1
            # insert
            length(new_codons) > mutator.max_codons && continue
            new_codon = random_codon(rng, gene_id_counter, geno.n_inputs)
            random_insert!(new_codons, new_codon)
        elseif choice == 2
            # delete
            length(new_codons) <= mutator.min_codons && continue
            index = rand(rng, 1:length(new_codons))
            deleteat!(new_codons, index)
        elseif choice == 3
            # move
            index = rand(rng, 1:length(new_codons))
            new_index = rand(rng, 1:length(new_codons))
            new_codons[index], new_codons[new_index] = new_codons[new_index], new_codons[index]
        end
    end
    # rewire 
    for _ in 1:mutator.num_swaps_per_mut
        index = rand(rng, 1:length(new_codons))
        new_id = count!(gene_id_counter)
        new_one, new_two = two_random_inputs(rng, geno.n_inputs)
        new_codons[index] = SortingNetworkCodon(new_id, new_one, new_two)
    end
    SortingNetworkGenotype(new_codons, geno.n_inputs)
end

function random_insert!(codons::Vector{SortingNetworkCodon}, codon::SortingNetworkCodon)
    index = rand(1:length(codons))
    insert!(codons, index, codon)
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
