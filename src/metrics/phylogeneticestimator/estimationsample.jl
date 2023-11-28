module PhylogeneticEstimator
using JLD2
using CoEvo
using CoEvo.States: State
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement
using ...Species.Phylogenetic: PhylogeneticSpecies
using ...Utils: save_statistical, display_stats, format_stat
using Statistics

using CoEvo.Metrics: Metric
Base.@kwdef struct PhylogeneticEstimatorMetric <: Metric
    name::String="PhylogeneticEstimatorMetric"
    path::String="data/archive.jld2"
    key::String="phylogeneticestimatorstats" # per-generation key
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
    @assert sum([typeof(s) <: PhylogeneticSpecies for s in state.species]) > 0 "No phylogenetic species in state"
    # Log metrics from all species
    measurement = PairwisePhylogeneticEstimationSampleMeasurement(
        Dict(s.id => s.measurements[PhylogeneticEstimationSampleMeasurement]
             for s in state.species
             if PhylogeneticEstimationSampleMeasurement in keys(s.measurements)
        )
    )
    measurement
end

function CoEvo.Archivers.archive!(
    archiver::CoEvo.Archivers.Basic.BasicArchiver, 
    report::CoEvo.Reporters.Basic.BasicReport{PhylogeneticEstimatorMetric,
                                              PairwisePhylogeneticEstimationSampleMeasurement}
)
    pairwise_measurements = report.measurement.measurements
    if report.to_print
        println("---")
        println("PhylogeneticEstimatorMetric")
        for speciesa_id in keys(pairwise_measurements)
            for speciesb_id in keys(pairwise_measurements[speciesa_id])
                m = pairwise_measurements[speciesa_id][speciesb_id]
                println("SpeciesA: $speciesa_id, SpeciesB: $speciesb_id")
                print("DistanceStatistics: ")
                display_stats(m.DistanceStatistics)
                println("ErrorStatistics")
                display_stats(m.ErrorStatistics)
                println("DistanceErrorCorrelation: $(format_stat(m.DistanceErrorCorrelation))")
                println("NumSamples: $(format_stat(m.NumSamples))")
            end
        end
    end
    if report.to_save
        met_path = isdir(archiver.archive_path) ? joinpath(archiver.archive_path, report.metric.path) : archiver.archive_path
        for speciesa_id in keys(pairwise_measurements)
            for speciesb_id in keys(pairwise_measurements[speciesa_id])
                m = pairwise_measurements[speciesa_id][speciesb_id]
                estimator_path = "gen/$(report.generation)/$(report.metric.key)/$speciesa_id-$speciesb_id"
                jldopen(met_path, "a+") do file
                    save_statistical(file, "$estimator_path/dist_stats", m.DistanceStatistics)
                    save_statistical(file, "$estimator_path/error_stats", m.ErrorStatistics)
                    file["$estimator_path/decorr"] = m.DistanceErrorCorrelation
                    file["$estimator_path/n_samples"] = m.NumSamples
                end
            end
        end
    end
end
end
