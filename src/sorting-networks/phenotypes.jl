export create_phenotype, act!, SortingNetworkPhenotypeCreator, SortingNetworkPhenotype, SortingNetworkTestCasePhenotype, SortingNetworkTestCasePhenotypeCreator, create_phenotype

using CoEvo.Phenotypes: Phenotype, PhenotypeCreator, act!, reset!

struct SortingNetworkPhenotypeCreator <: PhenotypeCreator
    n::Int64
end

struct SortingNetworkPhenotype <: Phenotype
    network::Array{Int64, 2}
    n::Int64
end

struct SortingNetworkTestCasePhenotype{N} <: Phenotype
    n::NTuple{N, Int64}
end
struct SortingNetworkTestCasePhenotypeCreator <: PhenotypeCreator
    n::Int64
end

function compare!(a::Int64, b::Int64, v::Array{Int64, 1})
    biggest, smallest = max(a, b), min(a, b)
    if v[biggest] < v[smallest]
        v[biggest], v[smallest] = v[smallest], v[biggest]
    end
end

netsort(snp::SortingNetworkPhenotype, sntc::SortingNetworkTestCasePhenotype) = netsort(snp, sntc.n)
function netsort(snp::SortingNetworkPhenotype, numbers::NTuple{N, Int64}) where {N}
    arr = [n for n in numbers]
    for (a, b) in eachrow(snp.network)
        compare!(a, b, arr)
    end
    return arr
end

function CoEvo.Phenotypes.create_phenotype(phenotype_creator::SortingNetworkPhenotypeCreator, geno::SortingNetworkGenotype)
    @assert phenotype_creator.n ∈ [2, 4, 8, 16]
    first_mask = UInt16(phenotype_creator.n - 1)  << 4
    second_mask = UInt16(phenotype_creator.n - 1)
    network = zeros(Int64, length(geno.codons), 2)
    num_active = 0
    for (i, codon) in enumerate(geno.codons)
        !is_active(codon) && continue
        num_active += 1
        network[num_active, 1] = 1+((codon.data & first_mask) >> 4) |> Int64
        network[num_active, 2] = 1+((codon.data & second_mask) >> 0) |> Int64
        @assert network[num_active, 1] ∈ 1:phenotype_creator.n
        @assert network[num_active, 2] ∈ 1:phenotype_creator.n
    end
    SortingNetworkPhenotype(network[1:num_active,:], phenotype_creator.n)
end

function CoEvo.Phenotypes.create_phenotype(phenotype_creator::SortingNetworkTestCasePhenotypeCreator, geno::SortingNetworkTestCaseGenotype)
    SortingNetworkTestCasePhenotype(geno.inputs)
end
