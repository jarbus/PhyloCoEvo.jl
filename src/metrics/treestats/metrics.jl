module TreeStatistics
using CoEvo
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement, GroupStatisticalMeasurement
using CoEvo.States: State
using CoEvo.Metrics: Metric
using ...Utils: save_statistical, display_stats, format_stat
using JLD2
using PhylogeneticTrees

Base.@kwdef struct TreeStatisticsMetric <: Metric
    name::String="TreeStatistics"
    path::String="data/archive.jld2"
    key::String="tree_stats" # per-generation key
end

struct IntGroupStatisticalMeasurement <: Measurement
    measurements::Dict{Int, BasicStatisticalMeasurement}
end

Base.@kwdef struct TreeStatisticsMeasurement <: Measurement
    n_nodes::Int
    per_distance_fitness_error_stats::IntGroupStatisticalMeasurement
    per_distance_interaction_error_stats::IntGroupStatisticalMeasurement
    total_distance_statistics::BasicStatisticalMeasurement
    distance_to_mrca::BasicStatisticalMeasurement
end

Base.@kwdef struct GroupTreeStatisticsMeasurement <: Measurement
    measurements::Dict{String, TreeStatisticsMeasurement}
end

function CoEvo.Genotypes.get_size(genotype::CoEvo.Genotypes.Vectors.BasicVectorGenotype)
    return length(genotype.genes)
end

function filter_pairwise_distances(pairwise_distances::Dict{Tuple{Int, Int}, Int})
    filtered_pairwise_distances = Dict{Tuple{Int, Int}, Int}()
    num_dists = zeros(Int, 11)
    for ((id1, id2), dist) in pairwise_distances
        dist > 5                   && continue
        num_dists[dist + 1] > 9999  && continue
        num_dists[dist + 1] += 1
        filtered_pairwise_distances[(id1, id2)] = dist
    end
    return filtered_pairwise_distances
end

function max_entries(dist_errors::Vector{Vector{Float64}})
    upper_limit = 9999
    dist_lens = length.(dist_errors)
    return all(dist_lens .> upper_limit)
end

function CoEvo.Metrics.measure(
    ::TreeStatisticsMetric,
    state::State
)
    # Assumes that each species interacts with only one other species at a time
    species_measurements = Dict{String, TreeStatisticsMeasurement}()

    species_pds = [filter_pairwise_distances(species.dist_data.pairwise_distances)
                   for species in state.species]
    dist_int_diffs = [Float64[] for _ in 0:10]
    for a in 1:(length(species_pds)-1), ((ind_a1, ind_a2), dist_a) in species_pds[a]
        max_entries(dist_int_diffs) && break
        
        ind_a1 ∉ keys(state.evaluations[1].outcomes) && continue
        ind_a2 ∉ keys(state.evaluations[1].outcomes) && continue
        outcomes_a1 = state.evaluations[1].outcomes[ind_a1]
        outcomes_a2 = state.evaluations[1].outcomes[ind_a2]
        for b in (a+1):length(species_pds), ((ind_b1,ind_b2), dist_b) in species_pds[b]
            dist = dist_a + dist_b
            dist > 10 && continue
            length(dist_int_diffs[dist+1]) > 9999 && continue
            # TODO clean this up
            if ind_b1 ∈ keys(outcomes_a1) && ind_b2 ∈ keys(outcomes_a2)
                outcome_1 = outcomes_a1[ind_b1]
                outcome_2 = outcomes_a2[ind_b2]
                estimation_error = abs(outcome_1 - outcome_2)
                push!(dist_int_diffs[dist+1], estimation_error)
            end

            if ind_b1 ∈ keys(outcomes_a2) && ind_b2 ∈ keys(outcomes_a1)
                outcome_1 = outcomes_a2[ind_b1]
                outcome_2 = outcomes_a1[ind_b2]
                estimation_error = abs(outcome_1 - outcome_2)
                push!(dist_int_diffs[dist+1], estimation_error)
            end
        end
    end
    
    for (species, evals) in zip(state.species, state.evaluations)
        mrca, pairwise_distances, mrca_distances = 
            species.dist_data.mrca,
            species.dist_data.pairwise_distances,
            species.dist_data.mrca_distances
        
        # Sample fitness differences for each distance
        dist_fit_diffs = [Float64[] for _ in 1:10]
        for rec1 in evals.records
            max_entries(dist_fit_diffs) && break
            for rec2 in evals.records
                id1, id2 = rec1.id, rec2.id
                id1 == id2 && continue
                (id1, id2) ∉ keys(pairwise_distances) && continue
                distance = pairwise_distances[id1, id2]
                distance > 10 && continue
                length(dist_fit_diffs[distance]) > 9999 && continue
                estimation_error = abs(rec1.fitness - rec2.fitness)
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
        for (distance, errors) in enumerate(dist_fit_diffs)
            length(errors) == 0 && continue
            per_distance_fitness_error_stats[distance] = BasicStatisticalMeasurement(errors)
        end
        # we include 0 interaction distance
        per_distance_interaction_error_stats = Dict{Int, BasicStatisticalMeasurement}()
        for (distance_plus_one, errors) in enumerate(dist_int_diffs)
            length(errors) == 0 && continue
            per_distance_interaction_error_stats[distance_plus_one-1] = BasicStatisticalMeasurement(errors)
        end
        total_distance_statistics = collect(values(pairwise_distances))
        mrca_distances = collect(values(mrca_distances))

        # Create TreeStatisticsMeasurement
        species_measurements[species.id] = TreeStatisticsMeasurement(
            n_nodes=n_nodes,
            per_distance_fitness_error_stats=IntGroupStatisticalMeasurement(
                per_distance_fitness_error_stats
            ),
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


function CoEvo.Archivers.archive!(
    archiver::CoEvo.Archivers.Basic.BasicArchiver, 
    report::CoEvo.Reporters.Basic.BasicReport{TreeStatisticsMetric, GroupTreeStatisticsMeasurement}
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
        met_path = isdir(archiver.archive_path) ? joinpath(archiver.archive_path, report.metric.path) : archiver.archive_path
        met_key = report.metric.key
        per_dist_int_err_stats = first(report.measurement.measurements).second.per_distance_interaction_error_stats.measurements
        jldopen(met_path, "a+") do file
            for (distance, stats) in per_dist_int_err_stats
                save_statistical(file, "gen/$(report.generation)/$met_key/dist_int_errors/$distance", stats)
            end
        end
    end
end
end
