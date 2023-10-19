module PhyloCoEvo
export PhylogeneticSpeciesCreator, PhylogeneticSpecies, create_species,
       EstimationPerformer, PhylogeneticMatchMaker, TreeStatisticsMetric,
       OutcomeScalarFitnessEvaluator, OutcomeScalarFitnessEvaluation,
       SortingNetworkGenotypeCreator, SortingNetworkGenotype
include("./species.jl")
include("./evaluation.jl")
include("./metrics.jl")
include("./estimation-performer.jl")
include("./matchmaker.jl")
include("./sort/networks/genotype/genotypes.jl")
end
