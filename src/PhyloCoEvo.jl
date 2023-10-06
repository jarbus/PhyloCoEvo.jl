module PhyloCoEvo
export PhylogeneticSpeciesCreator, PhylogeneticSpecies, create_species,
       DistanceError, EstimationPerformer, PhylogeneticMatchMaker, TreeStatisticsMetric
include("./phylogenetic.jl")
include("./metrics.jl")
include("./estimation-performer.jl")
include("./matchmaker.jl")
end
