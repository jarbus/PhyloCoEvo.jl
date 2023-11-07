using PhylogeneticTrees

@testset "Estimation" begin
    #   1
    #   |
    #   2
    #  / \
    # 3   4
    # |   |
    # 5   6
    #     |
    #     7
    tree = PhylogeneticTree([1])
    add_child!(tree, 1, 2)
    add_child!(tree, 2, 3)
    add_child!(tree, 2, 4)
    add_child!(tree, 3, 5)
    add_child!(tree, 4, 6)
    add_child!(tree, 6, 7)
    PhyloCoEvo.find_k_nearest_interactions
end
