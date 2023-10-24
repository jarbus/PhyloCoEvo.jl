export SortedMetric

Base.@kwdef struct SortedMetric <: SpeciesMetric
    name::String="Sorted"
    path::String="data/archive.jld2"
    key::String="sorted" # per-generation key
end

Base.@kwdef struct SortingNetworkMeasurement <: Measurement
    best_sorter::Float64 # percentage of 2^n inputs sorted correctly
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
        elseif genotype isa SortingNetworkTestCaseGenotype
            tc_fitnesses = collect(values(evaluation.fitnesses))
        else
            error("Unknown genotype $(typeof(species.pop[1].geno)) in SortedMetric")
        end
    end

    measurement = SortingNetworkMeasurement(
        sn_fitnesses=BasicStatisticalMeasurement(sn_fitnesses),
        tc_fitnesses=BasicStatisticalMeasurement(tc_fitnesses),
        best_sorter=best_score        
    )
    measurement
end

function CoEvo.archive!(
    archiver::BasicArchiver, 
    gen::Int, 
    report::BasicReport{SortedMetric, SortingNetworkMeasurement}
)
    # TODO finish this
    if report.to_print
        println("Archiving SortedMetric")
    end
    if report.to_save
    end
end
