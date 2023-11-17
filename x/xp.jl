#!/bin/env julia
using Revise
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true
ENV["JULIA_STACKTRACE_MINIMAL"] = true
using AbbreviatedStackTraces
using XPlot
using Plots

mp = "./23-10-23 sorting-network/make_plot.jl"
# list all files in src directory recursively
files = vcat(split( read(`find ../src -type f -name "*.jl"`, String), "\n"),
             split( read(`find ../test -type f -name "*.jl"`, String), "\n"), [mp, "./xp.jl"])
roc(files) do
    # cat make_plot.jl | julia
    include(mp)
end

