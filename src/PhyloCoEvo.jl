module PhyloCoEvo
export initialize_x
export PhylogeneticSpeciesCreator, PhylogeneticSpecies, create_species,
       EstimationPerformer, TreeStatisticsMetric,
       OutcomeScalarFitnessEvaluator, OutcomeScalarFitnessEvaluation,
       reset!
include("./utils.jl")
include("./species.jl")
include("./metrics.jl")
include("./estimators/estimators.jl")
using .Estimators: Estimator, Phylogenetic, estimate!
include("./evaluation.jl")
include("./estimation-performer.jl")
include("./matchmaker.jl")
include("./sorting-networks/evaluation.jl")
include("./sorting-networks/genotypes.jl")
include("./sorting-networks/phenotypes.jl")
include("./sorting-networks/domain.jl")
include("./sorting-networks/mutator.jl")
include("./sorting-networks/metrics.jl")
end
