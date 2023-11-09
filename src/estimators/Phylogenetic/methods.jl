function estimates_to_outcomes(estimates::Vector{EstimatedOutcome})
    ids = Set(id for e in estimates for id in (e.ida, e.idb))
    individual_outcomes = Dict{Int, Dict{Int, Float64}}(id=>Dict{Int, Float64}() for id in ids)
    for e in estimates
        individual_outcomes[e.ida][e.idb] = e.est_outcomea
        individual_outcomes[e.idb][e.ida] = e.est_outcomeb
    end
    sorted_dict_individual_outcomes = Dict(k=>SortedDict{Int,Float64}(v) for (k,v) in individual_outcomes)
    sorted_dict_individual_outcomes
end


function create_individual_outcomes_from_estimates(estimates::Vector{EstimatedOutcome})
    ids = Set(id for e in estimates for id in (e.ida, e.idb))
    individual_outcomes = Dict{Int, Dict{Int, Float64}}(id=>Dict{Int,Float64}() for id in ids)
    for e in estimates
        individual_outcomes[e.ida][e.idb] = e.est_outcomea
        individual_outcomes[e.idb][e.ida] = e.est_outcomeb
    end
    sorted_dict_individual_outcomes = Dict(k=>SortedDict{Int,Float64}(v) for (k,v) in individual_outcomes)
    sorted_dict_individual_outcomes
end

function measure_estimation_samples(estimates::Vector{EstimatedOutcome},
                         outcomes::Dict{Int, SortedDict{Int, Float64}})
    """Compute metrics of interest for a set of estimates"""
    distances = [d for e in estimates for d in e.distances]
    errorsa = [abs(e.est_outcomea - outcomes[e.ida][e.idb]) for e in estimates]
    errorsb = [abs(e.est_outcomeb - outcomes[e.idb][e.ida]) for e in estimates]
    PhylogeneticEstimationSampleMeasurement(distances, errorsa), PhylogeneticEstimationSampleMeasurement(distances, errorsb)
end

function two_layer_merge!(d1::Dict{Int, <:AbstractDict}, d2::Dict{Int, <:AbstractDict})
    """Merge dictionary of dictionaries `d2` into `d1` by merging the inner dictionaries
    if the key is in both dictionaries, and adding the key+dict if it is not in `d1`."""
    # TODO: Profile and Optimize
    for id in keys(d2)
        # If id is in individual_outcomes, merge the two dictionaries using merge!
        if id ∈ keys(d1)
            merge!(d1[id], d2[id])
        else
            @warn "Estimating all outcomes for individual $id"
            d1[id] = d2[id]
        end
    end
end

function weighted_average_outcome(related_outcomes::Vector{RelatedOutcome})
    """Compute the weighted average outcome of a set of related outcomes

    Arguments:
    =========
    related_outcomes: Vector{RelatedOutcome}

    Returns:
    ========
    weighted_average_a, weighted_average_b: Float64, Float64
        The weighted average outcome for individuals a and b
    """
    k = length(related_outcomes)
    if k == 1
        weights = [1.0]
    else
        dists = [related_outcomes[i].dist for i in 1:k]
        inv_dist =  sum(dists) .- dists
        weights = inv_dist ./ sum(inv_dist)
    end
    weighted_average_a, weighted_average_b = 0.0, 0.0
    for i in 1:k
        @inbounds weighted_average_a += related_outcomes[i].outcomea * weights[i]
        @inbounds weighted_average_b += related_outcomes[i].outcomeb * weights[i]
    end
    return weighted_average_a, weighted_average_b
end

