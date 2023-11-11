using PhyloCoEvo
using StableRNGs
using CoEvo
using Test
using CoEvo.Names
using CoEvo.Individuals.Basic: BasicIndividual
using CoEvo.Genotypes.Vectors: BasicVectorGenotype
using PhyloCoEvo.Species.Phylogenetic: PhylogeneticSpecies

function make_dummy_phylo_species(n_parents::Vector{Int}, n_children::Vector{Int}; first_gen::Bool=false)
    g = BasicVectorGenotype([0.0])
    parents, children, id = [], [], 1
    for (np, nc) in zip(n_parents, n_children)
        push!(parents, [BasicIndividual(i, g, Int[]) for i in id:id+np-1])
        id += np
        if first_gen
            # Children don't have parents for first gen
            push!(children, [BasicIndividual(i, g, Int[]) for i in id:id+nc-1])
        else
            push!(children, [BasicIndividual(i, g, Int[rand(id-np:id-1)]) for i in id:id+nc-1])
        end
        id += nc
    end
    species = [PhylogeneticSpecies("species$i", p, c) 
               for (i, (p,c)) in enumerate(zip(parents, children))]
end

@testset "PhyloCoEvo.jl" begin
    # Write your tests here.
    include("x/ng/1/config.jl")
    include("sorting-networks.jl")
    include("x/sn/1/config.jl")
    include("x/pg/1/config.jl")
    include("estimation.jl")
    include("matchmakers.jl")
    include("x/ng-mm/1/config.jl")
end
