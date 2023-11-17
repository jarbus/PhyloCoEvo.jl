module SortingNetwork
export SortedMetric

using CoEvo
using CoEvo.Metrics: Metric
using CoEvo.Measurements: Measurement
using CoEvo.Measurements.Statistical: BasicStatisticalMeasurement, GroupStatisticalMeasurement
using CoEvo.States: State
using CoEvo.Species: AbstractSpecies
using CoEvo.Phenotypes: create_phenotype
using LRUCache: LRU
using ...Evaluators.Outcome: OutcomeScalarFitnessEvaluation
using ...Genotypes.SortingNetwork: SortingNetworkGenotype, SortingNetworkTestCaseGenotype
using ...Phenotypes.SortingNetwork: SortingNetworkPhenotypeCreator, netsort
using ...Utils: save_statistical, display_stats, format_stat
using JLD2



Base.@kwdef struct SortedMetric <: Metric
    name::String="Sorted"
    path::String="data/archive.jld2"
    key::String="sorted" # per-generation key
    num_tests_per_parasite::Int
    hash_cache::LRU{Any, Float64} = LRU{Any,Float64}(maxsize=10000)
end

Base.@kwdef mutable struct SortingNetworkMeasurement <: Measurement
    """We define three categories of sorting networks:
    - best: the best sorting network found in the population
    - allpass: sorting networks that sort all parasites but not all possible inputs
    - perfect: sorting networks that sort all possible inputs
    """
    best_percent::Float64    # how many inputs are sorted by the best sorting network
    allpass_percent::Float64 # how many allpass sorting networks are there
    perfect_percent::Float64 # how many perfect sorting networks are there
    best_size::Int           # size of the best sorting network
    allpass_min_size::Int    # minimum size of allpass sorting networks
    perfect_min_size::Int    # minimum size of perfect sorting networks
    dist_allpass_sizes::BasicStatisticalMeasurement
    dist_perfect_sizes::BasicStatisticalMeasurement
    sn_fitnesses::BasicStatisticalMeasurement
    tc_fitnesses::BasicStatisticalMeasurement
end

function percent_sorted(genotype::SortingNetworkGenotype)
    snpc = SortingNetworkPhenotypeCreator(genotype.n_inputs)
    snp = CoEvo.Phenotypes.create_phenotype(snpc, genotype)
    n_sorted = 0
    n = 2^genotype.n_inputs
    for i in 1:n
        bits = [Int(d) for d in digits(i, base=2, pad=genotype.n_inputs)]
        sorted = sort(bits)
        if sorted == netsort(snp, bits)
            n_sorted += 1
        end
    end
    return n_sorted / n
end

function compute_network_metrics(metric::SortedMetric,
        species::AbstractSpecies,
        evaluation::OutcomeScalarFitnessEvaluation)

    # Compute best sorting network stats
    sn_fitnesses = [r.fitness for r in evaluation.records]
    best_id = evaluation.records[argmax(sn_fitnesses)].id
    pop_idx = findfirst(ind -> ind.id == best_id, species.population)
    if pop_idx == nothing 
        child_idx = findfirst(ind -> ind.id == best_id, species.children)
        best_genotype = species.children[child_idx].genotype
    else
        best_genotype = species.population[pop_idx].genotype
    end
    snpc = SortingNetworkPhenotypeCreator(best_genotype.n_inputs)
    network = x->create_phenotype(snpc, x).network

    best_score = get!(metric.hash_cache, network(best_genotype)) do
        percent_sorted(best_genotype)
    end
    best_comparators = length(best_genotype.codons)

    # Compute allpass sorting network stats
    sn_ids = [r.id for r in evaluation.records]
    allpass_ids = filter(sn_ids) do id
        all(values(evaluation.outcomes[id]) .> metric.num_tests_per_parasite)
    end
    allpass_inds = filter([species.population; species.children]) do ind
        ind.id in allpass_ids
    end
    allpass_sizes = [length(ind.genotype.codons) for ind in allpass_inds]
    allpass_min_size = length(allpass_sizes) > 0 ? minimum(allpass_sizes) : 0
    allpass_percent = length(allpass_inds) / length(sn_ids)
    allpass_scores = [get!(metric.hash_cache, network(ind.genotype)) do
        percent_sorted(ind.genotype)
    end for ind in allpass_inds]
    dist_allpass_sizes = BasicStatisticalMeasurement(allpass_sizes)
    
    # Compute perfect sorting network stats
    perfect_idxs = [score == 1.0 for score in allpass_scores]
    perfect_sizes = allpass_sizes[perfect_idxs]
    perfect_min_size = length(perfect_sizes) > 0 ? minimum(perfect_sizes) : 0
    perfect_percent = sum(perfect_idxs) / length(sn_ids)
    dist_perfect_sizes = BasicStatisticalMeasurement(perfect_sizes)

    return SortingNetworkMeasurement(
        best_score,
        allpass_percent,
        perfect_percent,
        best_comparators,
        allpass_min_size,
        perfect_min_size,
        dist_allpass_sizes,
        dist_perfect_sizes,
        BasicStatisticalMeasurement(sn_fitnesses),
        BasicStatisticalMeasurement(Float64[])
    )
end


function CoEvo.Metrics.measure(
    metric::SortedMetric,
    state::State
)
    measurement = nothing
    tc_fitnesses = Float64[]
    for (species, evaluation) in zip(state.species, state.evaluations)
        genotype = species.population[1].genotype
        if genotype isa SortingNetworkGenotype
            measurement = compute_network_metrics(metric, species, evaluation)

        elseif genotype isa SortingNetworkTestCaseGenotype
            tc_fitnesses = [r.fitness for r in evaluation.records]
        else
            error("Unknown genotype $(typeof(species.pop[1].geno)) in SortedMetric")
        end
    end
    @assert !isnothing(measurement)
    measurement.tc_fitnesses = BasicStatisticalMeasurement(tc_fitnesses)

    measurement
end

function CoEvo.Archivers.archive!(
    archiver::CoEvo.Archivers.Basic.BasicArchiver, 
    report::CoEvo.Reporters.Basic.BasicReport{SortedMetric, SortingNetworkMeasurement}
)
    gen=report.generation
    m = report.measurement

    if report.to_print
        println("----")
        println("SortedMetric:")
        for field in fieldnames(typeof(m))
            if getfield(m, field) isa BasicStatisticalMeasurement
                println("$(field):")
                display_stats(getfield(m, field))
            else
                println("$(field): $(getfield(m, field))")
            end
        end

    end
    if report.to_save
        met_path = isdir(archiver.archive_path) ? joinpath(archiver.archive_path, report.metric.path) : archiver.archive_path
        met_key = report.metric.key

        jldopen(met_path, "a+") do file
            for field in fieldnames(typeof(m))
                if getfield(m, field) isa BasicStatisticalMeasurement
                    save_statistical(file, "gen/$gen/$met_key/$field", getfield(m, field))
                else
                    file["gen/$gen/$met_key/$field"] = getfield(m, field)
                end
            end
        end
    end
end
end
