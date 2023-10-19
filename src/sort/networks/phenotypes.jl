export create_phenotype, act!

using CoEvo.Ecosystems.Species.Phenotypes.Abstract: Phenotype, PhenotypeCreator
using CoEvo.Ecosystems.Species.Phenotypes.Interfaces: act!, reset!

struct SortingNetworkPhenotypeCreator <: PhenotypeCreator
    n::Int64
end

struct SortingNetworkPhenotype <: Phenotype
    network::Array{Int64, 2}
    n::Int64
end

function CoEvo.create_phenotype(phenotype_creator::SortingNetworkPhenotypeCreator, geno::Genotype)::Phenotype
    # TODO: implement
    SortingNetworkPhenotype(zeros(Int64, phenotype_creator.n, phenotype_creator.n), phenotype_creator.n)
end

function CoEvo.act!(phenotype::Phenotype, ::Any)
end

function CoEvo.act!(phenotype::Phenotype)
    return act!(phenotype, nothing)
end

function CoEvo.Ecosystems.Species.Phenotypes.Interfaces.reset!(phenotype::Phenotype)
end
