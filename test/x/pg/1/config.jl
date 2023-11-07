using CoEvo
using CoEvo.Configurations.PredictionGame
using PhyloCoEvo

using CoEvo.Individuals: Basic.BasicIndividualCreator
using CoEvo.Replacers.Truncation: TruncationReplacer
using CoEvo.Selectors.FitnessProportionate: FitnessProportionateSelector
using CoEvo.Recombiners.Clone: CloneRecombiner
using CoEvo.Reporters
using CoEvo.Reporters.Runtime: RuntimeReporter
using CoEvo.Reporters.Basic: BasicReporter
using CoEvo.Selectors.Tournament: TournamentSelector
using CoEvo.Metrics.Genotypes: GenotypeSize
using CoEvo.Metrics.Evaluations: AllSpeciesFitness


XDIR = initialize_x(dirname(@__FILE__))
# Overwrite the default species
function CoEvo.Configurations.PredictionGame.make_species_creators(configuration::PredictionGameConfiguration)
    species_ids = make_species_ids(configuration)
    genotype_creator, phenotype_creator, mutators = make_substrate_types(configuration)
    evaluator, selector = make_reproducer_types(configuration)
    species_creators = [
        PhylogeneticSpeciesCreator(
            id = species_id,
            n_population = configuration.n_population,
            n_children = configuration.n_population,
            genotype_creator = genotype_creator,
            individual_creator = BasicIndividualCreator(),
            phenotype_creator = phenotype_creator,
            evaluator = evaluator,
            replacer = TruncationReplacer(n_truncate = configuration.n_population),
            selector = selector,
            recombiner = CloneRecombiner(),
            mutators = mutators,
        ) 
        for species_id in species_ids
    ]
    return species_creators
end

function CoEvo.Configurations.PredictionGame.make_archive_path(::PredictionGameConfiguration)
    return joinpath(XDIR, "data/archive.jld2")
end

function CoEvo.Configurations.PredictionGame.make_reporters(configuration::PredictionGameConfiguration)
    runtime_reporter = RuntimeReporter(print_interval = 1)
    reporters = Reporter[
        BasicReporter(metric = TreeStatisticsMetric(),
            save_interval = 1,
            print_interval = 1),
        # BasicReporter(
        #     metric = GenotypeSize(), 
        #     save_interval = save_interval, 
        #     print_interval = print_interval
        # ),
        # BasicReporter(
        #     metric = GenotypeSize(name = "MinimizedGenotypeSize", minimize = true),
        #     save_interval = save_interval,
        #     print_interval = print_interval
        # ),
        # BasicReporter(
        #     metric = AllSpeciesFitness(), 
        #     save_interval = save_interval, 
        #     print_interval = print_interval
        # ),
    ]
    return runtime_reporter, reporters
end

function make_reproducer_types(configuration::PredictionGameConfiguration)
    reproduction_method = configuration.reproduction_method
    if reproduction_method == :roulette
        evaluator = OutcomeScalarFitnessEvaluator()
        selector = FitnessProportionateSelector(n_parents = configuration.n_population)
    elseif reproduction_method == :disco
        evaluator = OutcomeNSGAIIEvaluator(
            maximize = true, perform_disco = true, max_clusters = configuration.max_clusters,
        )
        selector = TournamentSelector(
            n_parents = configuration.n_population, 
            tournament_size = configuration.tournament_size
        )
    else
        throw(ArgumentError("Unrecognized reproduction method: $reproduction_method"))
    end
    return evaluator, selector
end



configuration = PredictionGameConfiguration(
    substrate = :gnarl_networks,
    reproduction_method = :disco,
    game = :collision_game,
    ecosystem_topology = :three_species_mix,
    n_population = 100,
    communication_dimension = 1,
)

ecosystem = run!(configuration, n_generations = 2)
