using Random
using CoEvo.Individuals: Individual, Basic.BasicIndividualCreator
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

@testset "NumbersGameTest" begin

    XDIR = initialize_x(dirname(@__FILE__))
    function dummy_eco_creator(;
        id::String = "test",
        trial::Int = 1,
        rng::AbstractRNG = StableRNG(42),
        n_pop::Int = 2,
        species_id1::String = "a",
        species_id2::String = "b",
        interaction_id::String = "NumbersGame{Sum}",
        default_vector::Vector{Float64} = fill(0.0, 1),
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
                    genotype_creator = CoEvo.Genotypes.Vectors.BasicVectorGenotypeCreator(default_vector = default_vector),
                    phenotype_creator = DefaultPhenotypeCreator(),
                    individual_creator = BasicIndividualCreator(),
                    evaluator = OutcomeScalarFitnessEvaluator(),
                    replacer = GenerationalReplacer(n_elite = n_elite),
                    selector = FitnessProportionateSelector(n_parents = n_pop),
                    recombiner = CloneRecombiner(),
                    mutators = [BasicVectorMutator(noise_standard_deviation = 0.1)],
                ),
                PhylogeneticSpeciesCreator(
                    id = species_id2,
                    n_population = n_pop,
                    n_children = n_pop,
                    genotype_creator = CoEvo.Genotypes.Vectors.BasicVectorGenotypeCreator(default_vector = default_vector),
                    individual_creator = BasicIndividualCreator(),
                    phenotype_creator = DefaultPhenotypeCreator(),
                    evaluator = OutcomeScalarFitnessEvaluator(),
                    replacer = GenerationalReplacer(n_elite = n_elite),
                    selector = FitnessProportionateSelector(n_parents = n_pop),
                    recombiner = CloneRecombiner(),
                    mutators = [BasicVectorMutator(noise_standard_deviation = 0.1)],
                ),
            ],
            job_creator = BasicJobCreator(
                n_workers = 1,
                interactions = [
                    BasicInteraction(
                        id = interaction_id,
                        environment_creator = StatelessEnvironmentCreator(NumbersGameDomain(:Sum)),
                        species_ids = [species_id1, species_id2],
                        matchmaker = AllvsAllMatchMaker(),
                    ),
                ],
            ),
            performer = EstimationPerformer(n_workers = 1),
            individual_id_counter = BasicCounter(0),
            gene_id_counter = BasicCounter(0),
            state_creator = BasicCoevolutionaryStateCreator(),
            runtime_reporter = RuntimeReporter(),
            reporters = Reporter[
                # BasicReporter(metric = AllSpeciesFitness()),
                # BasicReporter(metric = GenotypeSum()),
                BasicReporter(metric = TreeStatisticsMetric(),
                              save_interval = 1,
                              print_interval = 1)
            ],
            archiver = BasicArchiver(archive_path = XDIR),
        )
        return eco_creator
    
    end
    
    eco_creator = dummy_eco_creator(n_pop = 100)
    eco = evolve!(eco_creator, n_generations=2)
    @test true
end
