using Random
using CoEvo.Individuals: Individual
using CoEvo.Genotypes.Vectors: BasicVectorGenotypeCreator
using CoEvo.Phenotypes.Defaults: DefaultPhenotypeCreator
using CoEvo.Replacers.Generational: GenerationalReplacer
using CoEvo.Selectors.FitnessProportionate: FitnessProportionateSelector
using CoEvo.Recombiners.Clone: CloneRecombiner
using CoEvo.Mutators.Vectors: BasicVectorMutator


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
            rng = rng,
            species_creators = [
                PhylogeneticSpeciesCreator(
                    id = species_id1,
                    n_pop = n_pop,
                    geno_creator = BasicVectorGenotypeCreator(default_vector = default_vector),
                    phenotype_creator = DefaultPhenotypeCreator(),
                    evaluator = OutcomeScalarFitnessEvaluator(),
                    replacer = GenerationalReplacer(n_elite = n_elite),
                    selector = FitnessProportionateSelector(n_parents = n_pop),
                    recombiner = CloneRecombiner(),
                    mutators = [BasicVectorMutator(noise_standard_deviation = 0.1)],
                ),
                PhylogeneticSpeciesCreator(
                    id = species_id2,
                    n_pop = n_pop,
                    geno_creator = BasicVectorGenotypeCreator(default_vector = default_vector),
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
                    interaction_id => BasicInteraction(
                        id = interaction_id,
                        environment_creator = StatelessEnvironmentCreator(NumbersGameDomain(:Sum)),
                        species_ids = [species_id1, species_id2],
                        matchmaker = PhylogeneticMatchMaker(type = :plus),
                    ),
                ],
            ),
            performer = EstimationPerformer(n_workers = 1),
            reporters = Reporter[
                # BasicReporter(metric = AllSpeciesFitness()),
                BasicReporter(metric = GenotypeSum()),
                BasicReporter(metric = TreeStatisticsMetric(),
                              save_interval = 1,
                              print_interval = 1)
            ],
            archiver = BasicArchiver(jld2_path = XDIR),
        )
        return eco_creator
    
    end
    
    eco_creator = dummy_eco_creator(n_pop = 100)
    eco = evolve!(eco_creator, n_gen=2)
    @test true
end
