using CoEvo.Measurements: BasicStatisticalMeasurement, GroupStatisticalMeasurement
Base.@kwdef struct DistanceError <: SpeciesMetric
    name::String="DistanceError"
end

function CoEvo.measure(
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


function CoEvo.archive!(
    ::BasicArchiver, 
    gen::Int, 
    report::BasicReport{DistanceError, GroupStatisticalMeasurement}
)
    for (species_id, measurement) in report.measurement.measurements
        println("----")
        println("Distance Error for species ", species_id)
        println("Error: ", measurement.mean)
    end
end
