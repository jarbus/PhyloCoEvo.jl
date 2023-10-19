@testset "SortingNetworks" begin
    sngc = SortingNetworkGenotypeCreator(1, 1)
    rng = StableRNG(1)
    gene_id_counter = Counter(0)
    n_pop = 1
    genotypes = CoEvo.create_genotypes(sngc, rng, gene_id_counter, n_pop)
end
