using Distributed
n_workers = 2
addprocs(n_workers-1)
using PProf
using Random
using PhyloCoEvo
using StableRNGs
using CoEvo.Individuals: Basic.BasicIndividualCreator
using CoEvo.Phenotypes.Defaults: DefaultPhenotypeCreator
using CoEvo.Replacers.Generational: GenerationalReplacer
using CoEvo.Selectors.FitnessProportionate: FitnessProportionateSelector
using CoEvo.Recombiners.Clone: CloneRecombiner
using CoEvo.Mutators.Vectors: BasicVectorMutator
using CoEvo.Domains.NumbersGame: NumbersGameDomain
using CoEvo.Environments.Stateless: StatelessEnvironmentCreator
using CoEvo.Interactions.Basic: BasicInteraction
using CoEvo.Jobs.Basic: BasicJobCreator
using CoEvo.Reporters: Reporter, Basic.BasicReporter, Runtime.RuntimeReporter
using CoEvo.Metrics.Genotypes: GenotypeSum
using CoEvo.Archivers.Basic: BasicArchiver
using CoEvo.Ecosystems: evolve!, Basic.BasicEcosystemCreator
using CoEvo.States.Basic: BasicCoevolutionaryStateCreator
using CoEvo.Counters.Basic: BasicCounter
using CoEvo.MatchMakers.AllvsAll: AllvsAllMatchMaker

XDIR = initialize_x(dirname(@__FILE__))
function dummy_eco_creator(;
    id::String = "test",
    trial::Int = 1,
    rng::AbstractRNG = StableRNG({SEED}),
    n_pop::Int = 2,
    species_id1::String = "a",
    species_id2::String = "b",
    interaction_id::String = "NumbersGame{Sum}",
    n_inputs::Int = 16,
    n_codons::Int = 80,
    n_elite::Int = 10
)
    eco_creator = BasicEcosystemCreator(
        id = id,
        trial = trial,
        random_number_generator = rng,
        species_creators = [
            PhylogeneticSpeciesCreator(
                id = species_id1,
                n_population = n_pop,
                n_children = n_pop,
                genotype_creator = SortingNetworkGenotypeCreator(n_codons, n_inputs),
                phenotype_creator = SortingNetworkPhenotypeCreator(n_inputs),
                individual_creator = BasicIndividualCreator(),
                evaluator = OutcomeScalarFitnessEvaluator(),
                replacer = GenerationalReplacer(n_elite = n_elite),
                selector = FitnessProportionateSelector(n_parents = n_pop),
                recombiner = CloneRecombiner(),
                mutators = [SortingNetworkMutator()],
            ),
            PhylogeneticSpeciesCreator(
                id = species_id2,
                n_population = 100,
                n_children = 100,
                genotype_creator = SortingNetworkTestCaseGenotypeCreator(20, n_inputs),
                phenotype_creator = SortingNetworkTestCasePhenotypeCreator(n_inputs),
                individual_creator = BasicIndividualCreator(),
                evaluator = OutcomeScalarFitnessEvaluator(),
                replacer = GenerationalReplacer(n_elite = n_elite),
                selector = FitnessProportionateSelector(n_parents = 100),
                recombiner = CloneRecombiner(),
                mutators = [SortingNetworkTestCaseMutator()],
            ),
        ],
        job_creator = BasicJobCreator(
            n_workers = n_workers,
            interactions = [
                BasicInteraction(
                    id = interaction_id,
                    environment_creator = StatelessEnvironmentCreator(SortingNetworkDomain(Partial())),
                    species_ids = [species_id1, species_id2],
                    matchmaker = AllvsAllMatchMaker(),
                ),
            ],
        ),
        performer = EstimationPerformer(n_workers = n_workers),
        individual_id_counter = BasicCounter(0),
        gene_id_counter = BasicCounter(0),
        state_creator = BasicCoevolutionaryStateCreator(),
        runtime_reporter = RuntimeReporter(),
        reporters = Reporter[
            # BasicReporter(metric = AllSpeciesFitness()),
            # BasicReporter(metric = GenotypeSum()),
            #BasicReporter(metric = TreeStatisticsMetric(),
            #              save_interval = 1,
            #              print_interval = 1)
            BasicReporter(metric=SortedMetric(),
                          save_interval=10,
                          print_interval=10)
        ],
        archiver = BasicArchiver(archive_path = XDIR),
    )
    return eco_creator

end

eco_creator = dummy_eco_creator(n_pop = 1000)
eco = evolve!(eco_creator, n_generations=2000)
