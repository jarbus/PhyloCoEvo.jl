using Revise
ENV["GKSwstype"]="nul"
using Glob
using XPlot
using Plots
default(size=(1200, 750), margin=30Plots.px)
mediadir = "./media"
# make mediadir if it doesnt exist
isdir(mediadir) || mkdir(mediadir)
function makeplot(files, met, savename; xlabel::String="Generation", ylabel::String=string(met), title=string(met))
    nc = NameConfig()
    # for plotting a separate plot for each file
    # xs = [agg(load(met, nc, pi)) for pi in files]
    # for plotting a single plot for all files
    xs = [agg(load(met, nc, files))]

    plots = [plot(x, ribbon=(0,0)) for x in xs]
    plot(plots..., xlabel=xlabel, ylabel=ylabel,title=title)
    figname=joinpath(mediadir, savename)
    savefig(figname)
    run(`kitten icat $figname`)
    plot()
end


makeplot(["."], EstimatorErrorStatistics() , "EstimatorErrorStatistics.png")
makeplot(["."], EstimatorDistanceStatistics(), "EstimatorDistanceStatistics.png")
makeplot(["."], DistanceErrorCorrelation(), "DistanceErrorCorrelation.png")

