using Revise
ENV["GKSwstype"]="nul"
using XPlot
using Plots
default(size=(1600, 900))
dir = "."
mediadir = "./media"
# make mediadir if it doesnt exist
isdir(mediadir) || mkdir(mediadir)
p = [
 ".",
]

nc = NameConfig()
iderrs = [agg(load(InteractionDistanceErrors(1:6), nc, pi)) for pi in p]
sfs = [agg(load(SortFits(), nc, pi)) for pi in p]
spc = [agg(load(BestSortPercentage(), nc, pi)) for pi in p]
ss = [agg(load(BestSortSize(), nc, pi)) for pi in p]

dist_plots = [plot(iderr, title=iderr[1].xname) for iderr in iderrs]
plot(dist_plots..., xlabel="Generation", ylabel="Error",)
figname=joinpath(mediadir, "distance_errors.png")
savefig(figname)
run(`kitten icat $figname`)
plot()

fitness_plots_network = plot(sfs[1][1], title=sfs[1][1].xname)
plot(fitness_plots_network, xlabel="Generation", ylabel="Fitness",)
figname=joinpath(mediadir, "fitnesses-net.png")
savefig(figname)
run(`kitten icat $figname`)
plot()


fitness_plots_tc = plot(sfs[1][2], title=sfs[1][2].xname)
plot(fitness_plots_tc, xlabel="Generation", ylabel="Fitness",)
figname=joinpath(mediadir, "fitnesses-tc.png")
savefig(figname)
run(`kitten icat $figname`)
plot()

percentage_plots = [plot(sp, title=sp[1].xname) for sp in spc]
plot(percentage_plots..., xlabel="Generation", ylabel="Percentage Solved",)
figname=joinpath(mediadir, "percentages.png")
savefig(figname)
run(`kitten icat $figname`)

size_plots = [plot(s, title=s[1].xname) for s in ss]
plot(size_plots..., xlabel="Generation", ylabel="Size",)
figname=joinpath(mediadir, "sizes.png")
savefig(figname)
run(`kitten icat $figname`)
