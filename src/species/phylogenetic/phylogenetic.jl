module Phylogenetic

using Random: AbstractRNG
using DataStructures: OrderedDict
using CoEvo
using PhylogeneticTrees
using CoEvo.Names
using CoEvo.Species: AbstractSpecies


struct PhylogeneticDistanceData
    mrca::Union{Int, Nothing}
    pairwise_distances::Dict{Tuple{Int,Int}, Int}
    mrca_distances::Dict{Int, Int}
end
PhylogeneticDistanceData() = PhylogeneticDistanceData(nothing, Dict(), Dict())
function PhylogeneticDistanceData(tree::PhylogeneticTree, ids::Set{Int}) 
    if rand() < 0.05 # purge unreachable nodes every twenty iterations or so
        mrca, pairwise_distances, mrca_distances = compute_pairwise_distances!(tree, ids; remove_unreachable_nodes=true, max_distance=11)
    else
        mrca, pairwise_distances, mrca_distances = nothing, Dict(), Dict()
    end
    PhylogeneticDistanceData(mrca, pairwise_distances, mrca_distances)
end


"""
    PhylogeneticSpecies{P <: PhenotypeCreator, I <: Individual}

Represents a species population and its offspring, with a phylogenetic tree.

# Fields
- `id::String`: Unique species identifier.
- `phenotype_creator::P`: Phenotype configuration.
- `pop::OrderedDict{Int, I}`: Current population.
- `children::OrderedDict{Int, I}`: Offspring of the population.
- `tree::PhylogeneticTree`: Phylogenetic tree of the population.
- `dist_data::PhylogeneticDistanceData`: Phylogenetic distance data.
- `randomly_sampled_interactions::Dict{String, Set{Tuple{Int,Int}}}`: Set of randomly sampled interactions between two species. Ids from this species are in the first position, and Ids from other species are in the second position
- `measurements::Dict{String,Any}`: measurements that are computed.
"""

struct PhylogeneticSpecies{I <: Individual} <: AbstractSpecies
    id::String
    population::Vector{I}
    children::Vector{I}
    tree::PhylogeneticTree
    dist_data::PhylogeneticDistanceData
    randomly_sampled_interactions::Dict{String, Set{Tuple{Int,Int}}}
    measurements::Dict{Any,Any}
end

function PhylogeneticSpecies(id::String,
        population::Vector{I},
        children::Vector{I},
        tree::PhylogeneticTree,
        dist_data::PhylogeneticDistanceData) where I <: Individual
    PhylogeneticSpecies(id, population, children, tree, dist_data, Dict{String,Set{Tuple{Int,Int}}}(), Dict{Any,Any}())
end

function PhylogeneticSpecies(id::String,
        population::Vector{I},
        children::Vector{I}) where I <: Individual
    parent_ids = [id for c in children for id in c.parent_ids]
    population_ids = [ind.id for ind in population]
    child_ids = [ind.id for ind in children]
    all_ids = [population_ids; child_ids]
    if length(parent_ids) == 0
        # If children have no parents (are genesis), make them genesis
        tree = PhylogeneticTree(all_ids)
    else
        # If children have parents, then add them to the tree as leaves
        tree = PhylogeneticTree(population_ids)
        for c in children
            add_child!(tree, c.parent_ids[1], c.id)
        end
    end

    dist_data = PhylogeneticDistanceData(tree, Set(all_ids))
    PhylogeneticSpecies(id, population, children, tree, dist_data)
end


end
