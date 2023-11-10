module Null
import ..Estimators: Estimator
struct NullEstimator <: Estimator end
using CoEvo.Species: AbstractSpecies
function estimate!(
    ::NullEstimator,
    ::Dict{Int, <:Dict{Int, Float64}},
    ::Vector{<:AbstractSpecies};)
end
end
