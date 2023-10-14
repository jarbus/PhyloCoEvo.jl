using CoEvo.Measurements: BasicStatisticalMeasurement, GroupStatisticalMeasurement
using JLD2
using PhylogeneticTrees

Base.@kwdef struct TreeStatisticsMetric <: SpeciesMetric
    name::String="TreeStatistics"
    path::String="data/TreeStatistics.jld2"
end

struct IntGroupStatisticalMeasurement <: Measurement
    measurements::Dict{Int, BasicStatisticalMeasurement}
end

Base.@kwdef struct TreeStatisticsMeasurement <: Measurement
    n_nodes::Int
    per_distance_fitness_error_stats::IntGroupStatisticalMeasurement
    per_distance_interaction_errors::Dict{Int, Vector{Float64}}
    per_distance_interaction_error_stats::IntGroupStatisticalMeasurement
    total_distance_statistics::BasicStatisticalMeasurement
    distance_to_mrca::BasicStatisticalMeasurement
end

Base.@kwdef struct GroupTreeStatisticsMeasurement <: Measurement
    measurements::Dict{String, TreeStatisticsMeasurement}
end

function filter_pairwise_distances(pairwise_distances::Dict{Tuple{Int, Int}, Int}, ids::Set{Int})
    filtered_pairwise_distances = Dict{Tuple{Int, Int}, Int}()
    for ((id1, id2), dist) in pairwise_distances
        (id1 <= id2) && dist <= 5 && id1 ∈ ids && id2 ∈ ids && (filtered_pairwise_distances[(id1, id2)] = pairwise_distances[(id1, id2)])
    end
    return filtered_pairwise_distances
end

function CoEvo.measure(
    ::Reporter{TreeStatisticsMetric},
    species_evaluations::Dict{<:AbstractSpecies, <:OutcomeScalarFitnessEvaluation},
    ::Vector{<:Observation}
)
    
    species_measurements = Dict{String, TreeStatisticsMeasurement}()

    @assert length(species_evaluations) == 2 "TreeStatisticsMetric only works for two species"


    dist_int_diffs = Dict{Int, Vector{Float64}}(i => Vector{Float64}() for i in 0:10)
    species = collect(keys(species_evaluations))
    species1_pop = species_evaluations[species[1]].fitnesses |> keys |> collect |> Set
    species2_pop = species_evaluations[species[2]].fitnesses |> keys |> collect |> Set
    species1_pd = filter_pairwise_distances(species[1].dist_data.pairwise_distances, species1_pop)
    species2_pd = filter_pairwise_distances(species[2].dist_data.pairwise_distances, species2_pop)
    for ((ind_a1,ind_a2), dist_a) in species1_pd
        outcomes_a1 = species_evaluations[species[1]].outcomes[ind_a1]
        outcomes_a2 = species_evaluations[species[1]].outcomes[ind_a2]
        for ((ind_b1,ind_b2), dist_b) in species2_pd
            dist = dist_a + dist_b
            outcome_1 = outcomes_a1[ind_b1]
            outcome_2 = outcomes_a2[ind_b2]
            estimation_error = abs(outcome_1 - outcome_2)
            push!(dist_int_diffs[dist], estimation_error)

            outcome_3 = outcomes_a1[ind_b2]
            outcome_4 = outcomes_a2[ind_b1]
            estimation_error = abs(outcome_3 - outcome_4)
            push!(dist_int_diffs[dist], estimation_error)
        end
    end
    # remove all empty distances
    for (dist, errors) in dist_int_diffs
        if length(errors) == 0
            delete!(dist_int_diffs, dist)
        end
    end
    
    for (species, evals) in species_evaluations
        mrca, pairwise_distances, mrca_distances = 
            species.dist_data.mrca,
            species.dist_data.pairwise_distances,
            species.dist_data.mrca_distances
        
        # Compute fitness differences for each distance
        dist_fit_diffs = Dict{Int, Vector{Float64}}()
        for (id1, fit1) in evals.fitnesses
            for (id2, fit2) in evals.fitnesses
                id1 == id2 && continue
                (id1, id2) ∉ keys(pairwise_distances) && continue
                distance = pairwise_distances[id1, id2]
                estimation_error = abs(fit1 - fit2)
                if !haskey(dist_fit_diffs, distance)
                    dist_fit_diffs[distance] = Vector{Float64}()
                end
                push!(dist_fit_diffs[distance], estimation_error)
            end
        end

        # Compute interaction differences for each distance
        # this is too costly
        # dist_int_diffs = Dict{Int, Vector{Float64}}()
        # for (inda1, outs1) in evals.outcomes, (indb1, out1) in outs1
        #     for (inda2, outs2) in evals.outcomes, (indb2, out2) in outs1
        #         (inda1, inda2) ∉ keys(pairwise_distances) && continue
        #         dist_across_species = pairwise_distances[inda1, inda2]
        #         found_other_distance = false
        #         for (other_species, other_evals) in species_evaluations
        #             other_species == species && continue
        #             @assert indb1 ∈ keys(other_evals.outcomes)
        #             @assert indb2 ∈ keys(other_evals.outcomes)
        #             (indb1, indb2) ∉ keys(other_species.dist_data.pairwise_distances) && break
        #             found_other_distance = true
        #             dist_across_species += other_species.dist_data.pairwise_distances[indb1, indb2]
        #         end
        #         estimation_error = abs(out1 - out2)
        #         if !haskey(dist_int_diffs, dist_across_species)
        #             dist_int_diffs[dist_across_species] = Vector{Float64}()
        #         end
        #         push!(dist_int_diffs[dist_across_species], estimation_error)
        #     end
        # end

        
        # Collect distributional data for TreeStatisticsMeasurement
        n_nodes = length(species.tree.tree)
        per_distance_fitness_error_stats = Dict{Int, BasicStatisticalMeasurement}()
        for (distance, errors) in dist_fit_diffs
            per_distance_fitness_error_stats[distance] = BasicStatisticalMeasurement(errors)
        end
        per_distance_interaction_error_stats = Dict{Int, BasicStatisticalMeasurement}()
        for (distance, errors) in dist_int_diffs
            per_distance_interaction_error_stats[distance] = BasicStatisticalMeasurement(errors)
        end
        total_distance_statistics = collect(values(pairwise_distances))
        mrca_distances = collect(values(mrca_distances))

        # Create TreeStatisticsMeasurement
        species_measurements[species.id] = TreeStatisticsMeasurement(
            n_nodes=n_nodes,
            per_distance_fitness_error_stats=IntGroupStatisticalMeasurement(
                per_distance_fitness_error_stats
            ),
            per_distance_interaction_errors=dist_int_diffs,
            per_distance_interaction_error_stats=IntGroupStatisticalMeasurement(
                per_distance_interaction_error_stats
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
    archiver::BasicArchiver, 
    gen::Int, 
    report::BasicReport{TreeStatisticsMetric, GroupTreeStatisticsMeasurement}
)
    if report.to_print
        for (species_id, measurement) in report.measurement.measurements
            println("----")
            println("Tree Statistics for species ", species_id)
            println("Number of nodes: ", measurement.n_nodes)

            print("Distance to MRCA:\n   ")
            display_stats(measurement.distance_to_mrca)

            println("Per-distance fitness statistics:")
            per_dist_fit_err_stats = measurement.per_distance_fitness_error_stats.measurements
            for distance in sort(collect(keys(per_dist_fit_err_stats)))
                print(format_stat(distance))
                display_stats(per_dist_fit_err_stats[distance])
            end

            println("Per-distance interaction statistics:")
            per_dist_int_err_stats = measurement.per_distance_interaction_error_stats.measurements
            for distance in sort(collect(keys(per_dist_int_err_stats)))
                print(format_stat(distance))
                display_stats(per_dist_int_err_stats[distance])
            end

            print("Total distance statistics:\n     ")
            display_stats(measurement.total_distance_statistics)
        end
    end
    if report.to_save
        # Log stats to jld2
        met_path = joinpath(archiver.jld2_path, report.metric.path)
        per_dist_int_err_stats = first(report.measurement.measurements).second.per_distance_interaction_error_stats.measurements
        # load existing data
        if gen == 1
            data = Dict("gen" => Dict())
        else
            data = load(met_path)
        end
        data["gen"]["$gen"] = Dict()
        for (distance, stats) in per_dist_int_err_stats
            data["gen"]["$gen"]["$distance"] = Dict()
            for field in fieldnames(typeof(stats))
                data["gen"]["$gen"]["$distance"]["$field"] = getfield(stats, field)
            end
        end
        # write data to tree-stats.jld2
        save(met_path, data)
        # Make violin plot of per_distance_interaction_errors
        # per_dist_int_err = first(report.measurement.measurements).second.per_distance_interaction_errors
        # if gen % 10 == 0
        #     plot_per_distance_errors(per_dist_int_err)
        # end
    end
end
