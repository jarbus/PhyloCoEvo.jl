module Phylogenetic

export estimate!

using ..Estimators: Estimator
using ...Species.Phylogenetic: PhylogeneticSpecies
using ...Metrics.PhylogeneticEstimator: PhylogeneticEstimationSampleMeasurement

import ..Estimators: estimate!    

using Statistics: mean, cor
using PhylogeneticTrees: PhylogeneticNode, PhylogeneticTree
using CoEvo.Species: AbstractSpecies
using DataStructures: SortedDict
using CoEvo.Observers: Observation
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement, GroupStatisticalMeasurement

include("./structs.jl")
include("./methods.jl")
include("./estimator.jl")
end
