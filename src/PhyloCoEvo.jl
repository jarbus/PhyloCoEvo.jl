module PhyloCoEvo
export initialize_x
export PhylogeneticSpeciesCreator, PhylogeneticSpecies, create_species,
       EstimationPerformer, TreeStatisticsMetric,
       OutcomeScalarFitnessEvaluator, OutcomeScalarFitnessEvaluation,
       reset!

using CoEvo

include("./utils.jl")
using .Utils

include("./species/species.jl")
using .Species: Species 

include("./species_creators/species_creators.jl")
using .SpeciesCreators: SpeciesCreators

include("./matchmakers/matchmakers.jl")
using .MatchMakers: MatchMakers

include("./evaluators/evaluators.jl")
using .Evaluators: Evaluators

include("./genotypes/genotypes.jl")
using .Genotypes: Genotypes

include("./phenotypes/phenotypes.jl")
using .Phenotypes: Phenotypes

include("./domains/domain.jl")
using .Domains: Domains

include("./mutators/mutator.jl")
using .Mutators: Mutators

include("./metrics/metrics.jl")
using .Metrics: Metrics

include("./selectors/selectors.jl")
using .Selectors: Selectors

include("./estimators/estimators.jl")
using .Estimators: Estimators

include("./ecosystems/ecosystems.jl")
# using .Ecosystems: Ecosystems

end
