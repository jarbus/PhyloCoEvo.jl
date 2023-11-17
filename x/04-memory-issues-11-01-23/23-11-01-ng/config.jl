n_workers=2
using Distributed
addprocs(n_workers-1)
using Random
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
using CoEvo.SpeciesCreators.Basic: BasicSpeciesCreator
using CoEvo.Performers.Basic: BasicPerformer
using CoEvo.Genotypes.Vectors: BasicVectorGenotypeCreator
using CoEvo.Replacers.Generational: GenerationalReplacer
using CoEvo.Evaluators.ScalarFitness: ScalarFitnessEvaluator
using CoEvo.Domains.NumbersGame: Focusing


id::String = "test"
trial::Int = 1
rng::AbstractRNG = StableRNG({SEED})
species_id1::String = "a"
species_id2::String = "b"
interaction_id::String = "NumbersGame{Focusing}"
default_vector::Vector{Float64} = [0.0]
n_elite=10
n_pop=1000


eco_creator = BasicEcosystemCreator(
    id = id,
    trial = trial,
    random_number_generator = rng,
    species_creators = [
             BasicSpeciesCreator(
                id = species_id1,
                n_population = n_pop,
                n_children = n_pop,
                genotype_creator = BasicVectorGenotypeCreator(default_vector = default_vector),
                phenotype_creator = DefaultPhenotypeCreator(),
                individual_creator = BasicIndividualCreator(),
                evaluator = ScalarFitnessEvaluator(),
                replacer = GenerationalReplacer(n_elite = n_elite),
                selector = FitnessProportionateSelector(n_parents = n_pop),
                recombiner = CloneRecombiner(),
                mutators = [BasicVectorMutator(noise_standard_deviation = 0.1)],
            ),
            BasicSpeciesCreator(
                id = species_id2,
                n_population = n_pop,
                n_children = n_pop,
                genotype_creator = BasicVectorGenotypeCreator(default_vector = default_vector),
                phenotype_creator = DefaultPhenotypeCreator(),
                individual_creator = BasicIndividualCreator(),
                evaluator = ScalarFitnessEvaluator(),
                replacer = GenerationalReplacer(n_elite = n_elite),
                selector = FitnessProportionateSelector(n_parents = n_pop),
                recombiner = CloneRecombiner(),
                mutators = [BasicVectorMutator(noise_standard_deviation = 0.1)],
            ),       
    ],
    job_creator = BasicJobCreator(
        n_workers = n_workers,
        interactions = [
            BasicInteraction(
                id = interaction_id,
                environment_creator = StatelessEnvironmentCreator(NumbersGameDomain(:Focusing)),
                species_ids = [species_id1, species_id2],
                matchmaker = AllvsAllMatchMaker(),
            ),
        ],
    ),
    performer = BasicPerformer(n_workers = n_workers),
    individual_id_counter = BasicCounter(0),
    gene_id_counter = BasicCounter(0),
    state_creator = BasicCoevolutionaryStateCreator(),
    runtime_reporter = RuntimeReporter(),
    reporters = Reporter[],
    archiver = BasicArchiver(archive_path = "."),
    garbage_collection_interval = 10,
)

eco = evolve!(eco_creator, n_generations=500)
