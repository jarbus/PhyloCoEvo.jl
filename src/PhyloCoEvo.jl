module PhyloCoEvo
export PhylogeneticSpeciesCreator, PhylogeneticSpecies, create_species,
       EstimationPerformer, PhylogeneticMatchMaker, TreeStatisticsMetric,
       OutcomeScalarFitnessEvaluator, OutcomeScalarFitnessEvaluation,
       reset!
include("./species.jl")
include("./evaluation.jl")
include("./metrics.jl")
include("./estimation-performer.jl")
include("./matchmaker.jl")
include("./sorting-networks/genotypes.jl")
include("./sorting-networks/phenotypes.jl")
include("./sorting-networks/domain.jl")
end
