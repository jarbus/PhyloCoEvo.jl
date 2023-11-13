export estimate!

function estimate!(
    estimator::Vector{Estimator},
    individual_outcomes::Dict{Int, <:AbstractDict{Int, Float64}},
    species::Vector{<:AbstractSpecies};)
    for estimator in estimator
        estimate!(estimator, individual_outcomes, species)
    end
end


function estimate!(
    estimator::Estimator,
    individual_outcomes::Dict{Int, <:AbstractDict{Int, Float64}},
    species::Vector{<:AbstractSpecies};)

    throw(ErrorException(
        "`estimate!` not implemented for $estimator and $species."))
end
