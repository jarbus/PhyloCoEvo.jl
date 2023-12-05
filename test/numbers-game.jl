using CoEvo.Names
using PhyloCoEvo.Domains.NumbersGame: CompareOnAll, CompareOnOne
using PhyloCoEvo.Mutators.Vectors: TwoIndexUniformVectorMutator
m = CoEvo.Domains.measure
@testset "CompareOnAll" begin
    domain = NumbersGameDomain(CompareOnAll())
    @test m(domain, [1., 2., 3.], [0., 1., 2.]) == [1., 0.]
    @test m(domain, [1., 1., 1.], [0., 1., 2.]) == [0., 0.]
    @test m(domain, [1., 1., 1.], [2., 2., 2.]) == [0., 1.]
    @test m(domain, [1., 1.], [1., 1.]) == [1., 1.]
    @test m(domain, [0., 4.], [1., 3.]) == [0., 0.]
end
@testset "CompareOnOne" begin
    domain = NumbersGameDomain(CompareOnOne())
    @test m(domain, [1., 2., 3.], [0., 1., 2.]) == [1., 0.]
    @test m(domain, [0., 1.5, 1.], [0., 1., 2.]) == [0., 0.]
    @test m(domain, [1., 1., 1.], [2., 2., 2.]) == [0., 1.]
    @test m(domain, [0., 4.], [1., 3.]) == [1., 0.]
    @test m(domain, [1., 1.], [1., 1.]) == [1., 1.]
end
@testset "TwoIndexUniformVectorMutator" begin
    mutator = TwoIndexUniformVectorMutator()
    rng = StableRNG(1)
    n_tests = 5
    for i in 1:n_tests
        len = rand(rng, 2:10)
        genotype = BasicVectorGenotype(randn(rng, len))
        mutated_genotype = mutate(mutator, rng, BasicCounter(), genotype) 
        # Confirm that only two indices were mutated
        @test sum(mutated_genotype.genes .== genotype.genes) == length(genotype) - 2
        # Confirm that the mutated genotype is within the bounds
        deltas = mutated_genotype.genes .- genotype.genes
        @test all(mutator.lower_bound .<= deltas .<= mutator.upper_bound)
    end
end
