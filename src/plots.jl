using StatsPlots
using Distributions

function plot_per_distance_errors(data::Dict{Int, Vector{Float64}})
    # each distance has a vector of errors
    xs = Vector{Int}()
    ys = Vector{Float64}()
    for (distance, errors) in data
        for error in errors
            push!(xs, distance)
            push!(ys, error)
        end
    end
    violin(xs, ys, legend = false)

    # plot normal distribution
    savefig("/tmp/violinplot.png")
    # display image in shell using icat kitten
    run(`kitten icat /tmp/violinplot.png`)
end
