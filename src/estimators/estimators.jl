module Estimators

export Estimator, Phylogenetic, estimate!
using CoEvo.Species: AbstractSpecies
include("abstract/abstract.jl")

include("interfaces/interfaces.jl")

include("phylogenetic/phylogenetic.jl")
using .Phylogenetic: Phylogenetic

include("null/null.jl")
using .Null: Null
end
