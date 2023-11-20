using StableRNGs
using PhyloCoEvo.Selectors.Lexicase: LexicaseSelector
using CoEvo.Selectors: select
using CoEvo.Individuals.Basic: BasicIndividual
using CoEvo.Evaluators.ScalarFitness: ScalarFitnessRecord
using PhyloCoEvo.Evaluators.Outcome: OutcomeScalarFitnessEvaluation
using CoEvo.Genotypes.Vectors: BasicVectorGenotype
using PhyloCoEvo.Selectors.Lexicase: fast_max_filter!

@testset "fast_max_filter!" begin
    source_ids = [1,2,3,4,5]
    n_source_ids = length(source_ids)
    target_ids = [0,0,0,0,0]
    outcomes = [0. 0. 0. 0. 0.;
                1. 1. 0. 0. 0.;
                1. 1. 1. 0. 0.;
                1. 1. 1. 1. 0.;
                1. 1. 1. 1. 1.]
    n_source_ids = fast_max_filter!(source_ids, n_source_ids, target_ids, outcomes, 1)
    @test target_ids[1:n_source_ids] == [2,3,4,5]
    source_ids, target_ids = target_ids, source_ids

    n_source_ids = fast_max_filter!(source_ids, n_source_ids, target_ids, outcomes, 4)
    @test target_ids[1:n_source_ids] == [4,5]
    source_ids, target_ids = target_ids, source_ids

    n_source_ids = fast_max_filter!(source_ids, n_source_ids, target_ids, outcomes, 3)
    @test target_ids[1:n_source_ids] == [4,5]
    source_ids, target_ids = target_ids, source_ids

    n_source_ids = fast_max_filter!(source_ids, n_source_ids, target_ids, outcomes, 5)
    @test target_ids[1:n_source_ids] == [5,]
end

@testset "LexicaseSelection" begin
    function get_parents(outcomes; n_parents=1000)
        s = LexicaseSelector(n_parents=100)
        rng = StableRNG(1)
        new_pop = [BasicIndividual(i, BasicVectorGenotype([0.0]), Int[]) for i in 1:5]
        evaluation = OutcomeScalarFitnessEvaluation("",ScalarFitnessRecord[], outcomes)
        [p.id for p in select(s, rng, new_pop, evaluation)]
    end
    @testset "All5s" begin
        # With this test, we expect 5 to be selected every single time
        outcomes = Dict(
            1 => Dict(6 => 1., 7 => 0., 8 => 0., 9 => 0., 10 => 0.),
            2 => Dict(6 => 1., 7 => 1., 8 => 0., 9 => 0., 10 => 0.),
            3 => Dict(6 => 1., 7 => 1., 8 => 1., 9 => 0., 10 => 0.),
            4 => Dict(6 => 1., 7 => 1., 8 => 1., 9 => 1., 10 => 0.),
            5 => Dict(6 => 1., 7 => 1., 8 => 1., 9 => 1., 10 => 1.)
        )
        @test all(get_parents(outcomes) .== 5)
    end
    @testset "5sAnd4s" begin
        # With this test, we expect 5 to be selected the most, followed by 4
        # followed by 3, 2 and 1 are never selected
        # If 9 or 10 are selected, 5 wins
        # If 8 is selected, 3 or 5 win
        # If 7 is selected, 4 wins
        # If 6 is selected, 3 or 4 win
        # So we expect:
        #   5 to be selected 50% of the time, 
        #   4 to be selected 30% of the time
        #   3 to be selected 20% of the time
        outcomes = Dict(
            1 => Dict(6 => 0., 7 => 0., 8 => 0., 9 => 0., 10 => 0.),
            2 => Dict(6 => 1., 7 => 0., 8 => 0., 9 => 0., 10 => 0.),
            3 => Dict(6 => 1., 7 => 0., 8 => 1., 9 => 0., 10 => 0.),
            4 => Dict(6 => 1., 7 => 1., 8 => 0., 9 => 0., 10 => 0.),
            5 => Dict(6 => 0., 7 => 0., 8 => 1., 9 => 1., 10 => 1.)
        )
        pids = get_parents(outcomes)
        @test sum(pids .== 5) > sum(pids .== 4) > sum(pids .== 3)
        @test sum(pids .== 2) == sum(pids .== 1) == 0
    end
    @testset "SpecialistAndGeneralist" begin
        # With this test, 1 is a specialist, 5 is a generalist, 
        # 2-4 are less fit generalists
        outcomes = Dict(
            1 => Dict(6 => 1., 7 => 0., 8 => 0., 9 => 0., 10 => 0.),
            2 => Dict(6 => 0., 7 => 1., 8 => 0.5, 9 => 1., 10 => 0.5),
            3 => Dict(6 => 0., 7 => 1., 8 => 0.5, 9 => 1., 10 => 0.5),
            4 => Dict(6 => 0., 7 => 1., 8 => 0.5, 9 => 1., 10 => 0.5),
            5 => Dict(6 => 0., 7 => 1., 8 => 1., 9 => 1., 10 => 1.)
        )
        pids = get_parents(outcomes)
        @test sum(pids .== 5) > sum(pids .== 1) 
        @test sum(pids .== 2) == sum(pids .== 3) == sum(pids .== 4) == 0
    end
end
