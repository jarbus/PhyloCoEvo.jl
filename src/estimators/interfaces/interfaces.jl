export estimate!

function estimate!(
    estimator::Vector{<:Estimator},
    individual_outcomes::Dict{Int, <:AbstractDict{Int, Float64}},
    species::Vector{<:AbstractSpecies};)
    for estimator in estimator
        estimate!(estimator, individual_outcomes, species)
    end

    # assert each individual_outcome has the same number of interactions
    unique_interaction_counts = unique(length.(values(individual_outcomes)))
    @assert length(unique_interaction_counts) == 1 "individual_outcomes have different number of interactions $(unique_interaction_counts)"
end


function estimate!(
    estimator::Estimator,
    individual_outcomes::Dict{Int, <:AbstractDict{Int, Float64}},
    species::Vector{<:AbstractSpecies};)

    throw(ErrorException(
        "`estimate!` not implemented for $estimator and $species."))
end
