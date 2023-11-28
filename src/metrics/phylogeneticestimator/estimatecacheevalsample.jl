module EstimateCacheEvalSample
using JLD2
using CoEvo
using CoEvo.States: State
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement
using ...Species.Phylogenetic: PhylogeneticSpecies
using ...Utils: save_statistical, display_stats, format_stat
using Statistics

using CoEvo.Metrics: Metric
Base.@kwdef struct EstimateCacheEvalSampleMetric <: Metric
    name::String="EstimateCacheEvalMetricSample"
    path::String="data/archive.jld2"
    key::String="estimatecacheevalsample" # per-generation key
end

Base.@kwdef struct EstimateCacheEvalSampleMeasurement <: Measurement
    NumEstimated::Int
    NumEvaluated::Int
    NumCached::Int
    NumSamples::Int
end

function CoEvo.Metrics.measure(
    ::EstimateCacheEvalSampleMetric,
    state::State
)
    # Aggregate the number of estimated, evaluated, cached, and sampled interactions
    # across all species in the state, and return a single measurement
    @assert sum([typeof(s) <: PhylogeneticSpecies for s in state.species]) > 0 "No phylogenetic species in state"
    total_num_estimated = 0
    total_num_evaluated = 0
    total_num_cached = 0
    total_num_samples = 0
    for s in state.species
        if EstimateCacheEvalSampleMeasurement in keys(s.measurements)
            for k in keys(s.measurements[EstimateCacheEvalSampleMeasurement])
                m = s.measurements[EstimateCacheEvalSampleMeasurement][k]
                total_num_estimated += m.NumEstimated
                total_num_evaluated += m.NumEvaluated
                total_num_cached += m.NumCached
                total_num_samples += m.NumSamples
            end
        end
    end
    # We divide by 2 because we are double-counting
    return EstimateCacheEvalSampleMeasurement(
        NumEstimated=total_num_estimated / 2,
        NumEvaluated=total_num_evaluated / 2,
        NumCached=total_num_cached / 2,
        NumSamples=total_num_samples / 2
    )
end

function CoEvo.Archivers.archive!(
    archiver::CoEvo.Archivers.Basic.BasicArchiver, 
    report::CoEvo.Reporters.Basic.BasicReport{EstimateCacheEvalSampleMetric, }
)
    m = report.measurement
    if report.to_print
        println("---")
        println("EstimateCacheEvaluationSample:")
        println("  NumEstimated: $(m.NumEstimated)")
        println("  NumEvaluated: $(m.NumEvaluated)")
        println("  NumCached: $(m.NumCached)")
        println("  NumSamples: $(m.NumSamples)")
    end
    if report.to_save
        met_path = isdir(archiver.archive_path) ? joinpath(archiver.archive_path, report.metric.path) : archiver.archive_path
        estimator_path = "gen/$(report.generation)/$(report.metric.key)"
        jldopen(met_path, "a+") do file
            file["$estimator_path/numestimated"] = m.NumEstimated
            file["$estimator_path/numevaluated"] = m.NumEvaluated
            file["$estimator_path/numcached"] = m.NumCached
            file["$estimator_path/numsamples"] = m.NumSamples
        end
    end
end

end
