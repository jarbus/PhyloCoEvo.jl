using Distributed
n_workers = 1
using Random
using PhyloCoEvo
using PhyloCoEvo.Evaluators.SortingNetwork: SortingNetworkEvaluator
using PhyloCoEvo.Evaluators.Outcome: OutcomeScalarFitnessEvaluator
using PhyloCoEvo.Phenotypes.SortingNetwork: SortingNetworkPhenotypeCreator
using PhyloCoEvo.Phenotypes.SortingNetwork: SortingNetworkTestCasePhenotypeCreator
using PhyloCoEvo.Genotypes.SortingNetwork: SortingNetworkGenotypeCreator
using PhyloCoEvo.Genotypes.SortingNetwork: SortingNetworkTestCaseGenotypeCreator
using PhyloCoEvo.Mutators.SortingNetwork: SortingNetworkMutator
using PhyloCoEvo.Mutators.SortingNetwork: SortingNetworkTestCaseMutator
using PhyloCoEvo.Domains.SortingNetwork: SortingNetworkDomain, Partial
using PhyloCoEvo.Metrics.SortingNetwork: SortedMetric
using PhyloCoEvo.SpeciesCreators.Phylogenetic: PhylogeneticSpeciesCreator
using CoEvo.Performers.Basic: BasicPerformer
using StableRNGs
using CoEvo.Individuals: Basic.BasicIndividualCreator
using CoEvo.Phenotypes.Defaults: DefaultPhenotypeCreator
using CoEvo.Replacers.Truncation: TruncationReplacer
using PhyloCoEvo.Selectors.Lexicase: LexicaseSelector
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
using PhyloCoEvo.MatchMakers: ParentsVsChildrenMatchMaker, RandomCohortMatchMaker
using PhyloCoEvo.Ecosystems.EstimatorEcosystem: EstimatorEcosystemCreator
using PhyloCoEvo.Estimators.Phylogenetic: PhylogeneticEstimator
using PhyloCoEvo.Estimators: Estimator
using PhyloCoEvo.Metrics.PhylogeneticEstimator: PhylogeneticEstimatorMetric

XDIR = initialize_x(dirname(@__FILE__))
id::String = "test"
trial::Int = 1
rng::AbstractRNG = StableRNG({SEED})
sorting_net_pop::Int = 500
sorting_net_parents::Int = 100
test_case_pop::Int = 500
test_case_parents::Int = 100
min_codons::Int = 60
max_codons::Int = 120
n_inputs::Int = 16
n_test_cases_per_parasite::Int = 10
species_id1::String = "a"
species_id2::String = "b"
interaction_id::String = "Sort"

eco_creator = EstimatorEcosystemCreator(basic= BasicEcosystemCreator(
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
            evaluator = SortingNetworkEvaluator(OutcomeScalarFitnessEvaluator(), n_test_cases_per_parasite),
            replacer = TruncationReplacer(n_truncate = sorting_net_pop),
            selector = LexicaseSelector(n_parents = sorting_net_parents,),
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
            replacer = TruncationReplacer(n_truncate = test_case_pop),
            selector = LexicaseSelector(n_parents = test_case_parents),
            recombiner = CloneRecombiner(),
            mutators = [SortingNetworkTestCaseMutator(num_swaps_per_mut = n_test_cases_per_parasite)],
        ),
    ],
    job_creator = BasicJobCreator(
        n_workers = n_workers,
        interactions = [
            BasicInteraction(
                id = interaction_id,
                environment_creator = StatelessEnvironmentCreator(SortingNetworkDomain(Partial())),
                species_ids = [species_id1, species_id2],
                matchmaker = RandomCohortMatchMaker(
                                n_matches_per_ind=Dict(
                                        species_id1=>test_case_parents,
                                        species_id2=>sorting_net_parents)
                                ,n_samples=500)
            ),
        ],
    ),
    performer = BasicPerformer(n_workers = n_workers),
    individual_id_counter = BasicCounter(0),
    gene_id_counter = BasicCounter(0),
    state_creator = BasicCoevolutionaryStateCreator(),
    runtime_reporter = RuntimeReporter(),
    reporters = Reporter[
        BasicReporter(metric=SortedMetric(num_tests_per_parasite=n_test_cases_per_parasite),
                      save_interval=10,
                      print_interval=10),
        BasicReporter(metric = PhylogeneticEstimatorMetric(),
                      save_interval = 5,
                      print_interval = 5),

    ],
    archiver = BasicArchiver(archive_path = joinpath(XDIR, "data/archive.jld2")),
    garbage_collection_interval = 10,
),
    estimators=[
        PhylogeneticEstimator(
            speciesa_id=species_id1,
            speciesb_id=species_id2,
            k=2,
            max_dist=10
        ),
    ],
)

eco = evolve!(eco_creator, n_generations=1000)
