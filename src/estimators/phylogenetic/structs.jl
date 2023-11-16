export PhylogeneticEstimator, PhylogeneticEstimationSampleMeasurement
"""
    struct PhylogeneticEstimator <: Estimator

The `PhylogeneticEstimator` struct represents a phylogenetic estimator used for coevolutionary analysis.

# Fields
- `k::Int`: The number of nearest interactions to use in estimation.
- `max_dist::Int`: The maximum interaction distance to look for before quiting.
"""
Base.@kwdef struct PhylogeneticEstimator <: Estimator
    speciesa_id::String
    speciesb_id::String
    k::Int
    max_dist::Int
    cached_outcomes::LRU{Int, Dict{Int, Float64}} = LRU{Int, Dict{Int, Float64}}(maxsize=100000)
end

function PhylogeneticEstimator(speciesa_id::String,
                               speciesb_id::String,
                               k::Int,
                               max_dist::Int)
    PhylogeneticEstimator(speciesa_id, speciesb_id, k, max_dist, LRU{Int, Dict{Int, Float64}}(maxsize=100000))
end


struct QueueElement
    indA::PhylogeneticNode
    indB::PhylogeneticNode
    dist::Int
    caller::Union{QueueElement, Nothing} # Prevent infinite loops by never returning to the caller
    add_B::Bool # Whether to add inds from tree B
end

struct RelatedOutcome
    """A related outcome is an interaction between two individuals, 
    and the distance to a query pair of individuals. Used for computing
    weighted averages of outcomes.
    """
    ida::Int
    idb::Int
    dist::Int
    outcomea::Float64
    outcomeb::Float64
end

struct EstimatedOutcome
    ida::Int
    idb::Int
    distances::Vector{Int} # the average distance of the k nearest interactions
    est_outcomea::Float64
    est_outcomeb::Float64
end

function EstimatedOutcome(ida::Int,
                          idb::Int,
                          related_outcomes::Vector{RelatedOutcome})
    @assert length(related_outcomes) > 0 "No related outcomes"
    wa, wb = weighted_average_outcome(related_outcomes)
    distances = [r.dist for r in related_outcomes]
    EstimatedOutcome(ida, idb, distances, wa, wb)
end



