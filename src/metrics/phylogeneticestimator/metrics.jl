module PhylogeneticEstimator
using CoEvo
using CoEvo.States: State
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement
using ...Species.Phylogenetic: PhylogeneticSpecies
using Statistics

using CoEvo.Metrics: Metric
Base.@kwdef struct PhylogeneticEstimatorMetric <: Metric
    name::String="TreeStatistics"
    path::String="data/archive.jld2"
    key::String="tree_stats" # per-generation key
end

Base.@kwdef struct PhylogeneticEstimationSampleMeasurement <: Measurement
    DistanceStatistics::BasicStatisticalMeasurement
    ErrorStatistics::BasicStatisticalMeasurement
    DistanceErrorCorrelation::Float64
    NumSamples::Int

    function PhylogeneticEstimationSampleMeasurement(distances::Vector{<:Real}, errors::Vector{<:Real})
        avg_distances = [mean(d) for d in distances]
        new(
            BasicStatisticalMeasurement(distances),
            BasicStatisticalMeasurement(errors),
            cor(avg_distances, errors),
            length(errors)
        )
    end
end

struct PairwisePhylogeneticEstimationSampleMeasurement <: Measurement
    measurements::Dict{String, Dict{String, PhylogeneticEstimationSampleMeasurement}}
end


function CoEvo.Metrics.measure(
    ::PhylogeneticEstimatorMetric,
    state::State
)
    @assert sum(s <: PhylogeneticSpecies for s in state.species) > 0 "No phylogenetic species in state"
    # Log metrics from all species
    PairwisePhylogeneticEstimationSampleMeasurement(
        Dict(s.id => s.measurements[PhylogeneticEstimationSampleMeasurement]
             for s in state.species)
       )
end

function CoEvo.Archivers.archive!(
    archiver::CoEvo.Archivers.Basic.BasicArchiver, 
    report::CoEvo.Reporters.Basic.BasicReport{PhylogeneticEstimatorMetric,
                                              PhylogeneticEstimationSampleMeasurement}
)
end
end
