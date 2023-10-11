using CoEvo.Measurements: BasicStatisticalMeasurement, GroupStatisticalMeasurement
using PhylogeneticTrees

Base.@kwdef struct TreeStatisticsMetric <: SpeciesMetric
    name::String="TreeStatistics"
end

struct IntGroupStatisticalMeasurement <: Measurement
    measurements::Dict{Int, BasicStatisticalMeasurement}
end

Base.@kwdef struct TreeStatisticsMeasurement <: Measurement
    n_nodes::Int
    per_distance_statistics::IntGroupStatisticalMeasurement
    total_distance_statistics::BasicStatisticalMeasurement
    distance_to_mrca::BasicStatisticalMeasurement
end

Base.@kwdef struct GroupTreeStatisticsMeasurement <: Measurement
    measurements::Dict{String, TreeStatisticsMeasurement}
end

function CoEvo.measure(
    ::Reporter{TreeStatisticsMetric},
    species_evaluations::Dict{<:AbstractSpecies, <:Evaluation},
    ::Vector{<:Observation}
)
    
    species_measurements = Dict{String, TreeStatisticsMeasurement}()

    for (species, evals) in species_evaluations
        ids = Set([collect(keys(species.pop));
                   collect(keys(species.children))])
        mrca, pairwise_distances, mrca_distances = 
            species.dist_data.mrca,
            species.dist_data.pairwise_distances,
            species.dist_data.mrca_distances
        
        distance_differences = Dict{Int, Vector{Float64}}()
        for (id1, fit1) in evals.fitnesses
            for (id2, fit2) in evals.fitnesses
                id1 == id2 && continue
                (id1, id2) ∉ keys(pairwise_distances) && continue
                distance = pairwise_distances[id1, id2]
                estimation_error = abs(fit1 - fit2)
                if !haskey(distance_differences, distance)
                    distance_differences[distance] = Vector{Float64}()
                end
                push!(distance_differences[distance], estimation_error)
            end
        end
        
        # Collect distributional data for TreeStatisticsMeasurement
        n_nodes = length(species.tree.tree)
        per_distance_statistics = Dict{Int, BasicStatisticalMeasurement}()
        for (distance, errors) in distance_differences
            per_distance_statistics[distance] = BasicStatisticalMeasurement(errors)
        end
        total_distance_statistics = collect(values(pairwise_distances))
        mrca_distances = collect(values(mrca_distances))

        # Create TreeStatisticsMeasurement
        species_measurements[species.id] = TreeStatisticsMeasurement(
            n_nodes=n_nodes,
            per_distance_statistics=IntGroupStatisticalMeasurement(
                per_distance_statistics
            ),
            total_distance_statistics=BasicStatisticalMeasurement(
                total_distance_statistics
            ),
            distance_to_mrca=BasicStatisticalMeasurement(
                mrca_distances
            ),
        )
    end
    measurement = GroupTreeStatisticsMeasurement(species_measurements)

    return measurement
end

format_stat(value::AbstractFloat) = lpad(string(round(value, digits=2)), 5)
format_stat(value::Int) = lpad(string(value), 5)

function display_stats(n::Int, min_value, mean_value, std_value, max_value)
    println("|$(format_stat(min_value))  $(format_stat(mean_value)) ± $(format_stat(std_value))  $(format_stat(max_value))| n=$n")
end

function display_stats(stat_measure::BasicStatisticalMeasurement)
    display_stats(
        stat_measure.n,
        stat_measure.minimum,
        stat_measure.mean,
        stat_measure.std,
        stat_measure.maximum,
    )
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

        print("Distance to MRCA:\n   ")
        display_stats(measurement.distance_to_mrca)

        println("Per-distance statistics:")
        per_dist_stats = measurement.per_distance_statistics.measurements
        for distance in sort(collect(keys(per_dist_stats)))
            print(format_stat(distance))
            display_stats(per_dist_stats[distance])
        end
        print("Total distance statistics:\n   ")
        display_stats(measurement.total_distance_statistics)
    end
end

