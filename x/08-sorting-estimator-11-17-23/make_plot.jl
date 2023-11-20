using Revise
ENV["GKSwstype"]="nul"
using XPlot
using Plots
default(size=(1200, 750), margin=30Plots.px)
dir = "."
mediadir = "./media"
# make mediadir if it doesnt exist
isdir(mediadir) || mkdir(mediadir)
p = [
 ".",
]

function makeplot(met, savename; xlabel::String="Generation", ylabel::String=string(met), title=string(met))
    nc = NameConfig()
    xs = [rolling(agg(load(met, nc, pi))) for pi in p]

    plots = [plot(x, ribbon=(0.0, 0.0)) for x in xs]
    plot(plots..., xlabel=xlabel, ylabel=ylabel,title=title)
    figname=joinpath(mediadir, savename)
    savefig(figname)
    run(`kitten icat $figname`)
    plot()
end

makeplot(BestSortSize(), "best_sort_size.png")
makeplot(AllPassPercentage(), "all_pass_percent.png")
makeplot(BestSortPercentage(), "best_sort_percent.png")
