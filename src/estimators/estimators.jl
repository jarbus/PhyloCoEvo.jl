module Estimators

export Estimator, Phylogenetic, estimate!
using CoEvo.Species: AbstractSpecies
include("abstract/abstract.jl")

include("interfaces/interfaces.jl")

include("Phylogenetic/phylogenetic.jl")
using .Phylogenetic: Phylogenetic
end
