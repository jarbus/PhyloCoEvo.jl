module PhyloCoEvo
export PhylogeneticSpeciesCreator, PhylogeneticSpecies, create_species,
       EstimationPerformer, PhylogeneticMatchMaker, TreeStatisticsMetric,
       OutcomeScalarFitnessEvaluator, OutcomeScalarFitnessEvaluation,
       plot_per_distance_errors
include("./species.jl")
include("./evaluation.jl")
include("./metrics.jl")
include("./estimation-performer.jl")
include("./matchmaker.jl")
include("./plots.jl")
end
