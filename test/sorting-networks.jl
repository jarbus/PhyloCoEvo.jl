@testset "SortingNetworks" begin
    sngc = SortingNetworkGenotypeCreator(1, 1, 1)
    snpc = SortingNetworkPhenotypeCreator(1)
    rng = StableRNG(1)
    gene_id_counter = Counter(0)
    n_pop = 1
    genotypes = create_genotypes(sngc, rng, gene_id_counter, n_pop)
    @test length(genotypes) == 1
    phenotype = create_phenotype(snpc, genotypes[1])
    domain = SortingNetworkDomain(SortingNetworkMetric())
end
