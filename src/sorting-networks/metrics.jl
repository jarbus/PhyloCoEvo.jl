export SortedMetric

Base.@kwdef struct SortedMetric <: Metric
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
    snp = PhyloCoEvo.Phenotypes.create_phenotype(snpc, genotype)
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


function CoEvo.Metrics.measure(
    ::SortedMetric,
    state::State
)
    best_score = 0
    best_comparators = 0
    best_n_inputs = 0
    sn_fitnesses = Float64[]
    tc_fitnesses = Float64[]
    for (species, evaluation) in zip(state.species, state.evaluations)
        genotype = species.population[1].genotype
        if genotype isa SortingNetworkGenotype
            sn_fitnesses = [r.fitness for r in evaluation.records]
            best_id = evaluation.records[argmax(sn_fitnesses)].id
            pop_idx = findfirst(ind -> ind.id == best_id, species.population)
            if pop_idx == nothing 
                child_idx = findfirst(ind -> ind.id == best_id, species.children)
                best_genotype = species.children[child_idx].genotype
            else
                best_genotype = species.population[pop_idx].genotype
            end

            best_score = percent_sorted(best_genotype)
            best_comparators = length(best_genotype.codons)
            best_n_inputs = best_genotype.n_inputs
        elseif genotype isa SortingNetworkTestCaseGenotype
            tc_fitnesses = [r.fitness for r in evaluation.records]
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

function CoEvo.Archivers.archive!(
    archiver::CoEvo.Archivers.Basic.BasicArchiver, 
    report::CoEvo.Reporters.Basic.BasicReport{SortedMetric, SortingNetworkMeasurement}
)
    gen=report.generation
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
        met_path = joinpath(archiver.archive_path, report.metric.path)
        met_key = report.metric.key
        sn_fit_stats = report.measurement.sn_fitnesses
        tc_fit_stats = report.measurement.tc_fitnesses

        jldopen(met_path, "a+") do file
            file["gen/$gen/$met_key/best_sorter_percent"] = report.measurement.best_sorter_percent
            file["gen/$gen/$met_key/best_sorter_size"] = report.measurement.best_sorter_size
            file["gen/$gen/$met_key/best_sorter_n_inputs"] = report.measurement.best_sorter_n_inputs
            save_statistical(file, "gen/$gen/$met_key/sn_fitnesses/", sn_fit_stats)
            save_statistical(file, "gen/$gen/$met_key/tc_fitnesses/", tc_fit_stats)
        end
    end
end
