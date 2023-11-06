using Distributed
n_workers = 1
using Random
using PhyloCoEvo
using StableRNGs
using CoEvo.Individuals: Basic.BasicIndividualCreator
using CoEvo.Phenotypes.Defaults: DefaultPhenotypeCreator
using CoEvo.Replacers.Truncation: TruncationReplacer
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
id::String = "test"
trial::Int = 1
rng::AbstractRNG = StableRNG({SEED})
sorting_net_pop::Int = 800
sorting_net_parents::Int = 10
sorting_net_elite::Int = 10
test_case_pop::Int = 200
test_case_parents::Int = 10
test_case_elite::Int = 10
min_codons::Int = 60
max_codons::Int = 120
n_inputs::Int = 16
n_test_cases_per_parasite::Int = 20
species_id1::String = "a"
species_id2::String = "b"
interaction_id::String = "Sort"

eco_creator = BasicEcosystemCreator(
    id = id,
    trial = trial,
    random_number_generator = rng,
    species_creators = [
        PhylogeneticSpeciesCreator(
            id = species_id1,
            n_population = sorting_net_pop,
            n_children = sorting_net_pop,
            genotype_creator = SortingNetworkGenotypeCreator(min_codons, n_inputs, min_codons,max_codons),
            phenotype_creator = SortingNetworkPhenotypeCreator(n_inputs),
            individual_creator = BasicIndividualCreator(),
            evaluator = OutcomeScalarFitnessEvaluator(),
            replacer = TruncationReplacer(n_truncate = sorting_net_elite),
            selector = FitnessProportionateSelector(n_parents = sorting_net_pop),
            recombiner = CloneRecombiner(),
            mutators = [SortingNetworkMutator()],
        ),
        PhylogeneticSpeciesCreator(
            id = species_id2,
            n_population = test_case_pop,
            n_children = test_case_pop,
            genotype_creator = SortingNetworkTestCaseGenotypeCreator(n_test_cases_per_parasite, n_inputs),
            phenotype_creator = SortingNetworkTestCasePhenotypeCreator(n_inputs),
            individual_creator = BasicIndividualCreator(),
            evaluator = OutcomeScalarFitnessEvaluator(),
            replacer = TruncationReplacer(n_truncate = test_case_elite),
            selector = FitnessProportionateSelector(n_parents = test_case_pop),
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
        BasicReporter(metric = TreeStatisticsMetric(),
                     save_interval = 5,
                     print_interval = 5)
        BasicReporter(metric=SortedMetric(),
                      save_interval=10,
                      print_interval=10)
    ],
    archiver = BasicArchiver(archive_path = XDIR),
    garbage_collection_interval = 10,
)

eco = evolve!(eco_creator, n_generations=1000)
