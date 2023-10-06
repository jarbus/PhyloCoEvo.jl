using CoEvo.Measurements: BasicStatisticalMeasurement, GroupStatisticalMeasurement
Base.@kwdef struct DistanceError <: SpeciesMetric
    name::String="DistanceError"
end

Base.@kwdef struct TreeStatisticsMetric <: SpeciesMetric
    name::String="TreeStatistics"
end

Base.@kwdef struct TreeStatisticsMeasurement <: Measurement
    n_nodes::Int
end

Base.@kwdef struct GroupTreeStatisticsMeasurement <: Measurement
    measurements::Dict{String, TreeStatisticsMeasurement}
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

function CoEvo.measure(
    ::Reporter{TreeStatisticsMetric},
    species_evaluations::Dict{<:AbstractSpecies, <:Evaluation},
    ::Vector{<:Observation}
)
    species_measurements = Dict(
        species.id => TreeStatisticsMeasurement(
            n_nodes=length(species.tree.tree)
        ) 
        for species in keys(species_evaluations)
    )
    measurement = GroupTreeStatisticsMeasurement(species_measurements)

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

function CoEvo.archive!(
    ::BasicArchiver, 
    gen::Int, 
    report::BasicReport{TreeStatisticsMetric, GroupTreeStatisticsMeasurement}
)
    for (species_id, measurement) in report.measurement.measurements
        println("----")
        println("Tree Statistics for species ", species_id)
        println("Number of nodes: ", measurement.n_nodes)
    end
end

