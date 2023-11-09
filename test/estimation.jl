using PhylogeneticTrees
using DataStructures: SortedDict
using PhyloCoEvo.Estimators.Phylogenetic: weighted_average_outcome, RelatedOutcome, find_k_nearest_interactions, two_layer_merge!, estimate!, compute_estimates, PhylogeneticEstimator
using CoEvo.Names

@testset "estimation.jl" begin
@testset "WeightedAverage" begin
    # Only the last three arguments matter for a RelatedOutcome:
    #   distance, outcome_a, outcome_b
    # two interactions which are equally close to the target interaction
    @assert weighted_average_outcome([RelatedOutcome(0,0,1,1,0),
                                      RelatedOutcome(0,0,1,0,1)]) == (0.5, 0.5)
    # Test one relative
    @assert weighted_average_outcome([RelatedOutcome(0,0,1,1,0)]) == (1, 0)

    # Test one close relative and one far relative
    @assert weighted_average_outcome([RelatedOutcome(0,0,1,1,0),
                                      RelatedOutcome(1,0,9,0,1)]) == (0.9, 0.1)
    # Test two close relatives and one far relative
    @assert weighted_average_outcome([RelatedOutcome(0,0,1,1,0),
                                      RelatedOutcome(0,0,1,1,0),
                                      RelatedOutcome(0,0,8,0,1)]) == (0.9, 0.1)
    # Test four far relatives
    @assert weighted_average_outcome([RelatedOutcome(0,0,8,1,0),
                                      RelatedOutcome(0,0,8,1,0),
                                      RelatedOutcome(0,0,8,1,0),
                                      RelatedOutcome(0,0,8,0,1)]) == (0.75, 0.25)
end
@testset "find_k_nearest_interactions" begin
    @testset "6,14" begin
    #   A             B
    #-------        ------
    #   1            10
    #   |             |
    #   2            11
    #  / \           / \
    # 3   4        12   13
    # |   |         |   |
    # 5   6        14   15
    treeA = PhylogeneticTree([1])
    add_child!(treeA, 1, 2)
    add_child!(treeA, 2, 3)
    add_child!(treeA, 2, 4)
    add_child!(treeA, 3, 5)
    add_child!(treeA, 4, 6)
    treeB = PhylogeneticTree([10])
    add_child!(treeB, 10, 11)
    add_child!(treeB, 11, 12)
    add_child!(treeB, 11, 13)
    add_child!(treeB, 12, 14)
    add_child!(treeB, 13, 15)
    
    # make empty interaction outcomes dict with all individuals in both trees
    io = Dict{Int, SortedDict{Int, Float64}}(i=>SortedDict{Int, Float64}() for i in 1:15)
    
    # first two interactions are close to (6,14)
    io[4][14],io[14][4] = 1, 0 
    io[6][12],io[12][6] = 1, 0
    # these two are far away
    io[1][10],io[10][1] = 0, 1
    io[2][11],io[11][2] = 0, 1

    expected_dists_from_6_14 = [1, 1, 4, 6]
    max_dist = maximum(expected_dists_from_6_14)
    for k in 1:4
        nearest = find_k_nearest_interactions(6, 14, treeA, treeB, io, k, max_dist=max_dist)
        @test length(nearest) == k
        @test [n.dist for n in nearest] == expected_dists_from_6_14[1:k]
    end
    # Test that it works when the trees are swapped
    for k in 1:4
        nearest = find_k_nearest_interactions(14, 6, treeB, treeA, io, k, max_dist=max_dist)
        @test length(nearest) == k
        @test [n.dist for n in nearest] == expected_dists_from_6_14[1:k]
    end

    end
    @testset "Disconnected" begin
    # Test that we don't find any interactions that are not reachable
    #  A            B
    #-----        -----
    # 1 2         3  4
    # | |         |  |
    # 5 6         7  8
    treeA = PhylogeneticTree([1, 2])
    add_child!(treeA, 1, 5)
    add_child!(treeA, 2, 6)
    treeB = PhylogeneticTree([3, 4])
    add_child!(treeB, 3, 7)
    add_child!(treeB, 4, 8)
    io = Dict{Int, SortedDict{Int, Float64}}(i=>SortedDict{Int, Float64}() for i in 1:8)

    io[1][3],io[3][1] = 1, 0
    io[2][4],io[4][2] = 1, 0
    io[2][8],io[8][2] = 1, 0

    k = 2
    nearest = find_k_nearest_interactions(5, 7, treeA, treeB, io, k, max_dist=5)
    @test length(nearest) == 1 # This raises a warning, but it's what we expect
    @test nearest[1] == RelatedOutcome(1, 3, 2, 1, 0)

    k = 3
    nearest = find_k_nearest_interactions(6, 8, treeA, treeB, io, k, max_dist=5)
    @test length(nearest) == 2 # This raises a warning, but it's what we expect
    @test nearest[1] == RelatedOutcome(2, 8, 1, 1, 0)
    @test nearest[2] == RelatedOutcome(2, 4, 2, 1, 0)
    end
end
    @testset "two_layer_merge!" begin
        d1 = Dict(0=>Dict(1=>0), 1=>Dict(2=>0, 3=>0))
        d2 = Dict(1=>Dict(4=>0), 2=>Dict(1=>0, 3=>0))
        # Raises a warning, but it's what we expect
        two_layer_merge!(d1, d2)
        @test d1 == Dict(0=>Dict(1=>0),
                         1=>Dict(2=>0, 3=>0, 4=>0),
                         2=>Dict(1=>0, 3=>0))
    end
    @testset "estimate!" begin
    # Make two species
    #  A             B
    #  1             3
    #  |             |
    #  2             4
    #
    #  Estimate 2v4 from: 1v4, 2v3, and 1v3
    #  1v4 : 1, 0, dist = 1 weight = 3 norm_weight = 3/8 = 0.375
    #  2v3 : 0, 1  dist = 1 weight = 3 norm_weight = 3/8 = 0.375
    #  1v3 : 0, 1  dist = 2 weight = 2 norm_weight = 2/8 = 0.25
    #  Expected weighted average outcome for 2v4: (0.375, 0.625)
    species = make_dummy_phylo_species([1, 1], [1, 1])
    new_individual_outcomes() = Dict(1=>SortedDict(4=>1.,3=>0.), 
             2=>SortedDict(3=>0.),
             3=>SortedDict(1=>1.,2=>1.),
             4=>SortedDict(1=>0.))
    phyloestimator = PhylogeneticEstimator(3, 10)
    individual_outcomes = new_individual_outcomes()
    estimate!(phyloestimator, individual_outcomes, species)
    @test 4 ∈ keys(individual_outcomes[2])
    @test 2 ∈ keys(individual_outcomes[4])
    @test individual_outcomes[2][4] == 0.375
    @test individual_outcomes[4][2] == 0.625

    phyloestimator = PhylogeneticEstimator(2, 10)
    individual_outcomes = new_individual_outcomes()
    estimate!(phyloestimator, individual_outcomes, species)
    @test 4 ∈ keys(individual_outcomes[2])
    @test 2 ∈ keys(individual_outcomes[4])
    @test individual_outcomes[2][4] == 0.5
    @test individual_outcomes[4][2] == 0.5
    end
end
