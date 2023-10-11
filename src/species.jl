using Random: AbstractRNG
using DataStructures: OrderedDict
using CoEvo
using PhylogeneticTrees


using CoEvo.Species.Abstract: AbstractSpecies, SpeciesCreator
using CoEvo.Species.Genotypes.Abstract: GenotypeCreator
using CoEvo.Species.Genotypes.Interfaces: create_genotypes
using CoEvo.Species.Phenotypes.Abstract: PhenotypeCreator
using CoEvo.Species.Evaluators.Abstract: Evaluator, Evaluation
using CoEvo.Species.Replacers.Abstract: Replacer
using CoEvo.Species.Replacers.Interfaces: replace
using CoEvo.Species.Selectors.Abstract: Selector
using CoEvo.Species.Selectors.Interfaces: select
using CoEvo.Species.Recombiners.Abstract: Recombiner
using CoEvo.Species.Recombiners.Interfaces: recombine
using CoEvo.Species.Mutators.Abstract: Mutator
using CoEvo.Species.Mutators.Interfaces: mutate
using CoEvo.Ecosystems.Utilities.Counters: next!


function add_children!(tree::PhylogeneticTree, children::Dict{Int, <:Individual})
    for (id, child) in children
        @assert id ∉ keys(tree.tree) "id: $id, keys: $(keys(tree.tree))"
        @assert length(child.parent_ids) == 1
        parent_id = child.parent_ids[1]
        add_child!(tree, parent_id, id)
    end
end

struct PhylogeneticDistanceData
    mrca::Union{Int, Nothing}
    pairwise_distances::Dict{Tuple{Int,Int}, Int}
    mrca_distances::Dict{Int, Int}
end
PhylogeneticDistanceData() = PhylogeneticDistanceData(nothing, Dict(), Dict())
function PhylogeneticDistanceData(tree::PhylogeneticTree, ids::Set{Int}) 
    mrca, pairwise_distances, mrca_distances = compute_pairwise_distances(tree, ids)
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
"""

struct PhylogeneticSpecies{I <: Individual} <: AbstractSpecies
    id::String
    pop::Dict{Int, I}
    children::Dict{Int, I}
    tree::PhylogeneticTree
    dist_data::PhylogeneticDistanceData
end

# Constructors
function PhylogeneticSpecies(
    id::String,
    pop::Vector{<:Individual},
    children::Vector{<:Individual},
    tree::PhylogeneticTree,
    dist_data::PhylogeneticDistanceData
)
    return PhylogeneticSpecies(
        id,
        Dict(indiv.id => indiv for indiv in pop),
        Dict(indiv.id => indiv for indiv in children),
        tree,
        dist_data
    )
end

function PhylogeneticSpecies(id::String, pop::Dict{Int, I}) where {I <: Individual}
    tree = PhylogeneticTree(collect(keys(pop)))
    dist_data = PhylogeneticDistanceData(tree, Set(collect(keys(pop))))
    return PhylogeneticSpecies(
        id, 
        pop, 
        Dict{Int, I}(), 
        tree,
        dist_data
    )
end


"""
    PhylogeneticSpeciesCreator{...}

Defines the parameters for species generation.

# Fields
- `id::String`: A unique identifier for the species.
- `n_pop::Int`: Size of the population.
- `geno_creator::G`: Genotype configuration.
- `phenotype_creator::P`: Phenotype configuration.
- `indiv_creator::I`: Individual configuration.
- `evaluator::E`: Evaluation configuration.
- `replacer::RP`: Mechanism for replacing old individuals with new ones.
- `selector::S`: Mechanism for selecting parents for reproduction.
- `recombiner::RC`: Mechanism for recombination (e.g., crossover).
- `mutators::Vector{M}`: A list of mutation mechanisms.
- `reporters::Vector{R}`: A list of reporters for gathering species metrics.
"""
Base.@kwdef struct PhylogeneticSpeciesCreator{
    G <: GenotypeCreator,
    P <: PhenotypeCreator,
    E <: Evaluator,
    RP <: Replacer,
    S <: Selector,
    RC <: Recombiner,
    M <: Mutator,
} <: SpeciesCreator
    id::String
    n_pop::Int
    geno_creator::G
    phenotype_creator::P
    evaluator::E
    replacer::RP
    selector::S
    recombiner::RC
    mutators::Vector{M}
end

"""
Generate a new population of individuals using genotype and phenotype configurations.
TODO: add pop to phylotree

# Arguments
- `creator::SpeciesCfg`: Creator for the species.
- `rng::AbstractRNG`: Random number generator.
- `indiv_id_counter::Counter`: Counter for generating unique individual IDs.
- `gene_id_counter::Counter`: Counter for generating unique gene IDs.
"""
function CoEvo.create_species(
    species_creator::PhylogeneticSpeciesCreator,
    rng::AbstractRNG, 
    indiv_id_counter::Counter = Counter(),
    gene_id_counter::Counter = Counter(),
)
    genos = create_genotypes(
        species_creator.geno_creator, rng, gene_id_counter, species_creator.n_pop
    ) 
    indiv_ids = next!(indiv_id_counter, species_creator.n_pop)
    pop = Dict(
        indiv_id => Individual(indiv_id, geno, Int[]) 
        for (indiv_id, geno) in zip(indiv_ids, genos)
    )
    return PhylogeneticSpecies(species_creator.id, pop)
end

"""
Core reproduction phase of the evolutionary algorithm.

# Arguments
- `creator::SpeciesCfg`: Creator for the species.
- `rng::AbstractRNG`: Random number generator.
- `indiv_id_counter::Counter`: Counter for generating unique individual IDs.
- `gene_id_counter::Counter`: Counter for generating unique gene IDs.
- `species::Species`: Current species.
- `results::Vector{<:InteractionResult`: Interaction results of the individuals.

# Returns
- A new `PhylogeneticSpecies` containing the next generation population and their children.
"""
function CoEvo.create_species(
    species_creator::PhylogeneticSpeciesCreator,
    rng::AbstractRNG, 
    indiv_id_counter::Counter,  
    gene_id_counter::Counter,  
    species::PhylogeneticSpecies,
    evaluation::Evaluation
) 
    new_pop = replace(species_creator.replacer, rng, species, evaluation)
    parents = select(species_creator.selector, rng, new_pop, evaluation)
    new_children = recombine(species_creator.recombiner, rng, indiv_id_counter, parents)
    for mutator in species_creator.mutators
        new_children = mutate(mutator, rng, gene_id_counter, new_children)
    end
    new_children = Dict(indiv.id => indiv for indiv in new_children)
    # Update tree and compute distance data
    add_children!(species.tree, new_children)
    ids = Set([collect(keys(new_pop)); collect(keys(new_children))])
    dist_data = PhylogeneticDistanceData(species.tree, ids)

    new_species = PhylogeneticSpecies(species_creator.id, new_pop, new_children, species.tree, dist_data)
    return new_species
end
