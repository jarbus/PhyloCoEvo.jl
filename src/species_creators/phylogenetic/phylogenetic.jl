module Phylogenetic

using Random
using CoEvo
using CoEvo.Names
using CoEvo.Evaluators: Evaluation
using ..Species.Phylogenetic: PhylogeneticSpecies
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
    G <: CoEvo.Genotypes.GenotypeCreator,
    I <: CoEvo.Individuals.IndividualCreator,
    P <: CoEvo.Phenotypes.PhenotypeCreator,
    E <: CoEvo.Evaluators.Evaluator,
    RP <: CoEvo.Replacers.Replacer,
    S <: CoEvo.Selectors.Selector,
    RC <: CoEvo.Recombiners.Recombiner,
    M <: CoEvo.Mutators.Mutator,
} <: CoEvo.SpeciesCreators.SpeciesCreator
    id::String
    n_population::Int
    n_children::Int
    genotype_creator::G
    individual_creator::I
    phenotype_creator::P
    evaluator::E
    replacer::RP
    selector::S
    recombiner::RC
    mutators::Vector{M}
end

function CoEvo.SpeciesCreators.create_species(
    species_creator::PhylogeneticSpeciesCreator,
    random_number_generator::AbstractRNG, 
    individual_id_counter::Counter,
    gene_id_counter::Counter
)
    population = create_individuals(
        species_creator.individual_creator, 
        random_number_generator, 
        species_creator.genotype_creator, 
        species_creator.n_population, 
        individual_id_counter, 
        gene_id_counter
    )
    children = create_individuals(
        species_creator.individual_creator, 
        random_number_generator, 
        species_creator.genotype_creator, 
        species_creator.n_children, 
        individual_id_counter, 
        gene_id_counter
    )

    PhylogeneticSpecies(species_creator.id, population, children)
end

function CoEvo.SpeciesCreators.create_species(
    species_creator::PhylogeneticSpeciesCreator,
    random_number_generator::AbstractRNG, 
    individual_id_counter::Counter,  
    gene_id_counter::Counter,  
    species::PhylogeneticSpecies,
    evaluation::Evaluation
) 
    new_population = replace(
        species_creator.replacer, random_number_generator, species, evaluation
    )
    parents = select(
        species_creator.selector, random_number_generator, new_population, evaluation
    )
    new_children = recombine(
        species_creator.recombiner, random_number_generator, individual_id_counter, parents
    )
    for mutator in species_creator.mutators
        new_children = mutate(mutator, random_number_generator, gene_id_counter, new_children)
    end

    # Update tree and compute distance data
    add_children!(species.tree, new_children)
    ids = Set(vcat([ind.id for ind in new_population],[child.id for child in new_children]))
    dist_data = PhylogeneticDistanceData(species.tree, ids)
    
    new_species = PhylogeneticSpecies(species_creator.id,new_population, new_children, species.tree, dist_data)
    return new_species
end


function CoEvo.Ecosystems.Basic.construct_new_species(
    state::BasicCoevolutionaryState, species_creators::Vector{<:PhylogeneticSpeciesCreator}
)
    new_species = [
        CoEvo.SpeciesCreators.create_species(
            species_creators[index],
            state.random_number_generator, 
            state.individual_id_counter,
            state.gene_id_counter,
            state.species[index],
            state.evaluations[index]
        ) for (index) in eachindex(species_creators)
    ]
    return new_species
end
end
