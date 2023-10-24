export initialize_x
using CoEvo: BasicStatisticalMeasurement
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
