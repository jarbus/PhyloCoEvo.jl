module Metrics

include("./treestats/metrics.jl")
using .TreeStatistics

include("./sorting-network/metrics.jl")
using .SortingNetwork

include("./phylogeneticestimator/metrics.jl")
using .PhylogeneticEstimator

end
