export initialize_x
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
