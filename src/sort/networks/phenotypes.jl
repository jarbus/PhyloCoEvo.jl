export create_phenotype, act!, SortingNetworkPhenotypeCreator, SortingNetworkPhenotype, create_phenotype

using CoEvo.Ecosystems.Species.Phenotypes.Abstract: Phenotype, PhenotypeCreator
using CoEvo.Ecosystems.Species.Phenotypes.Interfaces: act!, reset!

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

netsort(numbers::NTuple{N, Int64}, snp::SortingNetworkPhenotype) where {N} = netsort(snp, numbers)
function netsort(snp::SortingNetworkPhenotype, numbers::NTuple{N, Int64}) where {N}
    arr = [n for n in numbers]
    for (a, b) in eachrow(snp.network)
        compare!(a, b, arr)
    end
    return arr
end

function CoEvo.create_phenotype(phenotype_creator::SortingNetworkPhenotypeCreator, geno::Genotype)::Phenotype
    # TODO: implement
    SortingNetworkPhenotype(zeros(Int64, phenotype_creator.n, phenotype_creator.n), phenotype_creator.n)
end

function CoEvo.act!(snp::SortingNetworkPhenotype, tcp::SortingNetworkTestCasePhenotype)
end

function CoEvo.act!(phenotype::Phenotype)
    # return act!(phenotype, nothing)
end

function CoEvo.Ecosystems.Species.Phenotypes.Interfaces.reset!(phenotype::SortingNetworkPhenotype)
end
