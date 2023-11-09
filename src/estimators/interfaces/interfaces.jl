export estimate!

function estimate!(
    estimator::Estimator,
    individual_outcomes::Dict{Int, <:Dict{Int, Float64}},
    species::Vector{<:AbstractSpecies};)

    throw(ErrorException(
        "`estimate!` not implemented for $estimator and $species."))
end
