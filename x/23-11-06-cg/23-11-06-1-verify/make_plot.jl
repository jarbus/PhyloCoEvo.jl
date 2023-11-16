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

nc = NameConfig()
iderrs = [rolling(agg(load(InteractionDistanceErrors(1:6), nc, pi))) for pi in p]

dist_plots = [plot(iderr, title=iderr[1].xname, ribbon=(0.0, 0.0)) for iderr in iderrs]
plot(dist_plots..., xlabel="Generation", ylabel="Error",title="Rolling average of distance errors")
figname=joinpath(mediadir, "distance_errors.png")
savefig(figname)
run(`kitten icat $figname`)
plot()

