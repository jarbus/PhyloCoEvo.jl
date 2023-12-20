using Revise
ENV["GKSwstype"]="nul"
using Glob
using XPlot
using Plots
default(size=(600, 400), margin=10Plots.px)
mediadir = "./media"
# make mediadir if it doesnt exist
isdir(mediadir) || mkdir(mediadir)

function merge_eval(evaluated, metric)
    data = XPlot.AggregatedTimeSeriesDataPoint[]
    num_evals = 0
    idx = 1
    for point in evaluated.data
        num_evals += point.value
        metpoint = metric.data[idx]
        if point.x == metpoint.x
            new_point = XPlot.AggregatedTimeSeriesDataPoint(num_evals, metpoint.value, metpoint.lower_bound, metpoint.upper_bound, metpoint.count)
            push!(data, new_point)
            idx += 1
        end
    end
    XPlot.AggregatedTimeSeriesData(metric.name, data, metric.xname, metric.yaxis, metric.label)
end

function makeplot(files, met, savename=string(met); ylabel::String=string(met), title=string(met))
    nc = NameConfig()
    # for plotting a separate plot for each file
    # xs = [agg(load(met, nc, pi)) for pi in files]
    xs = agg(load(Evaluated(), nc, files))
    # for plotting a single plot for all files
    ys = agg(load(met, nc, files))

    p = plot([merge_eval(x,y) for (x,y) in zip(xs,ys)])
    plot(p, xlabel="Number Evaluations", ylabel=ylabel,title=title)
    figname=joinpath(mediadir, savename)
    savefig(figname * ".png")
    savefig(figname * ".pdf")
    run(`kitten icat $figname.png`)
    plot()
end


makeplot(["."], PerfectPercentage())
# makeplot(["."], BestSortPercentage())
makeplot(["."], BestSortSize())
# makeplot(glob("0[1-2]*"), EstimatorErrorStatistics())

