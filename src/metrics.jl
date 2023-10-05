using CoEvo: GroupStatisticalMeasurement
Base.@kwdef struct DistanceError <: SpeciesMetric
    name::String="DistanceError"
end

function measure(
    ::Reporter{DistanceError},
    species_evaluations::Dict{<:AbstractSpecies, <:Evaluation},
    ::Vector{<:Observation}
)
    species_measurements = Dict(
        species.id => BasicStatisticalMeasurement(
            [sum(individual.geno.genes) for individual in values(species.pop)]
        ) 
        for species in keys(species_evaluations)
    )
    measurement = GroupStatisticalMeasurement(species_measurements)
    return measurement
end
