using Revise
ENV["GKSwstype"]="nul"
using XPlot
using Plots
dir = "."
mediadir = "./media"
# make mediadir if it doesnt exist
isdir(mediadir) || mkdir(mediadir)
p = [
 ".",
]

nc = NameConfig()
iderrs = [agg(load(InteractionDistanceErrors(1:5), nc, pi)) for pi in p]
sfs = [agg(load(SortFits(), nc, pi)) for pi in p]
spc = [agg(load(BestSortPercentage(), nc, pi)) for pi in p]

dist_plots = [plot(iderr, title=iderr[1].xname) for iderr in iderrs]
plot(dist_plots..., xlabel="Generation", ylabel="Error",)
figname=joinpath(mediadir, "distance_errors.png")
savefig(figname)
run(`kitten icat $figname`)

fitness_plots = [plot(sf, title=sf[1].xname) for sf in sfs]
plot(fitness_plots..., xlabel="Generation", ylabel="Fitness",)
figname=joinpath(mediadir, "fitnesses.png")
savefig(figname)
run(`kitten icat $figname`)

percentage_plots = [plot(sp, title=sp[1].xname) for sp in spc]
plot(percentage_plots..., xlabel="Generation", ylabel="Percentage Solved",)
figname=joinpath(mediadir, "percentages.png")
savefig(figname)
run(`kitten icat $figname`)
