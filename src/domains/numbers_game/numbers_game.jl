module NumbersGame
"""
The Compare-on-all and Compare-on-one Metrics for the Numbers Game Domain, from the original DISCO paper:
Liskowski, Paweł, and Krzysztof Krawiec. Online Discovery of Search Objectives for Test-Based Problems. Sept. 2017.
"""

export CompareOnAll, CompareOnOne
using CoEvo
using CoEvo.Domains.NumbersGame: NumbersGameDomain

import CoEvo.Domains: measure
using CoEvo.Metrics: Metric
using CoEvo.Domains: Domain
using Base: @kwdef


@kwdef struct CompareOnAll <: Metric 
    name::String = "CompareOnAll"
end
@kwdef struct CompareOnOne <: Metric 
    name::String = "CompareOnOne"
end

function CoEvo.Domains.measure(::NumbersGameDomain{CompareOnAll}, A::Vector{<:Real}, B::Vector{<:Real})
    # Check if all elements of A are greater than all elements of B, or vice versa
    # If so, return [1.0, 0.0] or [0.0, 1.0], respectively, otherwise return [0.0, 0.0]
    Float64[all(A .>= B), all(A .<= B)]
end

function CoEvo.Domains.measure(::NumbersGameDomain{CompareOnOne}, A::Vector{<:Real}, B::Vector{<:Real})
    # Gives A a score of 1 if it's element at B's max dimension is >= B's max value, and vice versa
    max_a_dim = argmax(A) 
    max_b_dim = argmax(B)
    a_score = A[max_b_dim] >= B[max_b_dim] ? 1.0 : 0.0
    b_score = B[max_a_dim] >= A[max_a_dim] ? 1.0 : 0.0
    return [a_score, b_score]

end

end
