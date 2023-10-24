export SortedMetric

Base.@kwdef struct SortedMetric <: SpeciesMetric
    name::String="Sorted"
    path::String="data/archive.jld2"
    key::String="sorted" # per-generation key
end

Base.@kwdef struct SortingNetworkMeasurement <: Measurement
    best_sorter_percent::Float64 # percentage of 2^n inputs sorted correctly
    best_sorter_size::Int # number of comparators in the best sorting network
    best_sorter_n_inputs::Int # number of inputs to the best sorting network
    sn_fitnesses::BasicStatisticalMeasurement
    tc_fitnesses::BasicStatisticalMeasurement
end

function percent_sorted(genotype::SortingNetworkGenotype)
    snpc = SortingNetworkPhenotypeCreator(genotype.n_inputs)
    snp = create_phenotype(snpc, genotype)
    n_sorted = 0
    n = 2^genotype.n_inputs
    for i in 1:n
        bits = [Int(d) for d in digits(i, base=2, pad=genotype.n_inputs)]
        sorted = sort(bits)
        if sorted == netsort(snp, Tuple(bits))
            n_sorted += 1
        end
    end
    return n_sorted / n
end


function CoEvo.measure(
    ::Reporter{SortedMetric},
    species_evaluations::Dict{<:AbstractSpecies, <:OutcomeScalarFitnessEvaluation},
    ::Vector{<:Observation}
)
    best_score = 0
    best_comparators = 0
    best_n_inputs = 0
    sn_fitnesses = Float64[]
    tc_fitnesses = Float64[]
    for (species, evaluation) in species_evaluations
        genotype = first(species.pop)[2].geno
        if genotype isa SortingNetworkGenotype
            sn_fitnesses = collect(values(evaluation.fitnesses))
            best_fitness = maximum(sn_fitnesses)
            best_id = findfirst(v->v==best_fitness, evaluation.fitnesses)
            best_genotype = best_id ∈ keys(species.pop) ? species.pop[best_id].geno : species.children[best_id].geno
            best_score = percent_sorted(best_genotype)
            best_comparators = num_active(best_genotype)
            best_n_inputs = best_genotype.n_inputs
        elseif genotype isa SortingNetworkTestCaseGenotype
            tc_fitnesses = collect(values(evaluation.fitnesses))
        else
            error("Unknown genotype $(typeof(species.pop[1].geno)) in SortedMetric")
        end
    end

    measurement = SortingNetworkMeasurement(
        sn_fitnesses=BasicStatisticalMeasurement(sn_fitnesses),
        tc_fitnesses=BasicStatisticalMeasurement(tc_fitnesses),
        best_sorter_percent=best_score,     
        best_sorter_size=best_comparators,
        best_sorter_n_inputs=best_n_inputs
    )
    measurement
end

function CoEvo.archive!(
    archiver::BasicArchiver, 
    gen::Int, 
    report::BasicReport{SortedMetric, SortingNetworkMeasurement}
)
    if report.to_print
        println("----")
        println("SortedMetric:")
        println("Best sorter percentage: $(report.measurement.best_sorter_percent)")
        println("Best sorter size: $(report.measurement.best_sorter_size)")
        print("Sorted fitnesses:\n  ")
        display_stats(report.measurement.sn_fitnesses)
        print("Test case fitnesses:\n  ")
        display_stats(report.measurement.tc_fitnesses)
    end
    if report.to_save
        met_path = joinpath(archiver.jld2_path, report.metric.path)
        met_key = report.metric.key
        sn_fit_stats = report.measurement.sn_fitnesses
        tc_fit_stats = report.measurement.tc_fitnesses

        jldopen(met_path, "a+") do file
            file["gen/$gen/$met_key/best_sorter_percent"] = report.measurement.best_sorter_percent
            file["gen/$gen/$met_key/best_sorter_size"] = report.measurement.best_sorter_size
            file["gen/$gen/$met_key/best_sorter_n_inputs"] = report.measurement.best_sorter_n_inputs
            save_statistical(file, "gen/$gen/$met_key/sn_fitnesses/$field", sn_fit_stats)
            save_statistical(file, "gen/$gen/$met_key/tc_fitnesses/$field", tc_fit_stats)
        end
    end
end
