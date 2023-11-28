module Phylogenetic

export estimate!

using ..Estimators: Estimator
using ...Species.Phylogenetic: PhylogeneticSpecies
using ...Metrics.PhylogeneticEstimator: PhylogeneticEstimationSampleMeasurement
using ...Metrics.EstimateCacheEvalSample: EstimateCacheEvalSampleMeasurement

import ..Estimators: estimate!    

using Statistics: mean, cor
using PhylogeneticTrees: PhylogeneticNode, PhylogeneticTree
using CoEvo.Species: AbstractSpecies
using DataStructures: SortedDict
using CoEvo.Observers: Observation
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement, GroupStatisticalMeasurement
using LRUCache: LRU

include("./structs.jl")
include("./methods.jl")
include("./estimator.jl")
end
