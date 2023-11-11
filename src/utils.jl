module Utils
export initialize_x, save_statistical, two_rand
using Random
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement
function initialize_x(path::String)
    """Takes in an experiment directory as input, creates data directory"""
    @assert isdir(path) "Path $(path) does not exist"
    if !isdir(joinpath(path, "data"))
        mkdir(joinpath(path, "data"))
    else 
        rm(joinpath(path, "data"), recursive=true)
        mkdir(joinpath(path, "data"))
    end
    path
end

function save_statistical(file, path::String,stats::BasicStatisticalMeasurement)
    for field in fieldnames(typeof(stats))
        file[joinpath(path, String(field))] = getfield(stats, field)
    end
end

function two_rand(rng::AbstractRNG, r::UnitRange)
    """Returns two random numbers from range r"""
    @assert r.start < r.stop "Range $(r) is not valid"
    a = rand(rng, r)
    b = rand(rng, r)
    while a == b
        b = rand(rng, r)
    end
    return a, b
end

format_stat(value::AbstractFloat) = lpad(string(round(value, digits=2)), 5)
format_stat(value::Int) = lpad(string(value), 5)

function display_stats(n::Int, min_value, mean_value, std_value, max_value)
    println("|$(format_stat(min_value))  $(format_stat(mean_value)) ± $(format_stat(std_value))  $(format_stat(max_value))| n=$n")
end

function display_stats(stat_measure::BasicStatisticalMeasurement)
    display_stats(
        stat_measure.n_samples,
        stat_measure.minimum,
        stat_measure.mean,
        stat_measure.std,
        stat_measure.maximum,
    )
end
end
