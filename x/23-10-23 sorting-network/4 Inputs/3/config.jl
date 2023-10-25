using Random
using CoEvo
using PhyloCoEvo
using StableRNGs

XDIR = initialize_x(dirname(@__FILE__))
function dummy_eco_creator(;
    id::String = "test",
    trial::Int = 1,
    rng::AbstractRNG = StableRNG(11121),
    n_pop::Int = 2,
    species_id1::String = "a",
    species_id2::String = "b",
    interaction_id::String = "SortingNetworks",
    n_elite::Int = 10
)
    eco_creator = BasicEcosystemCreator(
        id = id,
        trial = trial,
        rng = rng,
        species_creators = Dict(
            species_id1 => PhylogeneticSpeciesCreator(
                id = species_id1,
                n_pop = n_pop,
                geno_creator = SortingNetworkGenotypeCreator(6, 4),
                phenotype_creator = SortingNetworkPhenotypeCreator(4),
                evaluator = OutcomeScalarFitnessEvaluator(),
                replacer = GenerationalReplacer(n_elite = n_elite),
                selector = FitnessProportionateSelector(n_parents = n_pop),
                recombiner = CloneRecombiner(),
                mutators = [SortingNetworkMutator()],
            ),
            species_id2 => PhylogeneticSpeciesCreator(
                id = species_id2,
                n_pop = n_pop,
                geno_creator = SortingNetworkTestCaseGenotypeCreator(4),
                phenotype_creator = SortingNetworkTestCasePhenotypeCreator(4),
                evaluator = OutcomeScalarFitnessEvaluator(),
                replacer = GenerationalReplacer(n_elite = n_elite),
                selector = FitnessProportionateSelector(n_parents = n_pop),
                recombiner = CloneRecombiner(),
                mutators = [SortingNetworkTestCaseMutator()],
            ),
        ),
        job_creator = BasicJobCreator(
            n_workers = 1,
            interactions = Dict(
                interaction_id => BasicInteraction(
                    id = interaction_id,
                    environment_creator = StatelessEnvironmentCreator(SortingNetworkDomain(Partial())),
                    species_ids = [species_id1, species_id2],
                    matchmaker = PhylogeneticMatchMaker(type = :plus),
                ),
            ),
        ),
        performer = EstimationPerformer(n_workers = 1),
        reporters = Reporter[
            BasicReporter(metric = AllSpeciesFitness()),
            # BasicReporter(metric = GenotypeSum()),
            BasicReporter(metric=SortedMetric(), save_interval=1, print_interval=1),
            BasicReporter(metric = TreeStatisticsMetric(),
                          save_interval = 1,
                          print_interval = 1)
        ],
        archiver = BasicArchiver(jld2_path = XDIR),
    )
    return eco_creator

end

eco_creator = dummy_eco_creator(n_pop = 100)
eco = evolve!(eco_creator, n_gen=100)
