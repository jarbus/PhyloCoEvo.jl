using CoEvo.Individuals.Basic: BasicIndividual
using CoEvo.Genotypes.Vectors: BasicVectorGenotype
using PhyloCoEvo.MatchMakers: RandomCohortMatchMaker, ParentsVsAllMatchMaker
function test_random_cohort_matches(pop_sizes::Vector{Int},
                                    n_matches_per_ind::Vector{Int};
                                    n_samples::Int=0,
                                    gen::Int=2,
                                    n_gens_before_sampling::Int=1,
                                    n_expected_samples::Int=0,
                                    )
    "Function goes over each species pair and tests that the correct number of matches are made.
    The number of matches made between two species is the same for all species pairs.
    The number of matches made by an individual is the same for all individuals in a species.
    "
    @assert length(pop_sizes) == length(n_matches_per_ind)
    # We can handle two species with different pop sizes, or n species with the same pop size
    @assert length(pop_sizes) == 2 || length(Set(pop_sizes)) == 1
    # With random cohort matches, the number of matches made between two species is the same
    # for all species pairs.
    @assert length(Set(ps * nm for (ps, nm) in zip(pop_sizes, n_matches_per_ind))) == 1
    # each pair has the same num of matches
    n_species_pairs = length(pop_sizes) * (length(pop_sizes) - 1) / 2 
    if gen > n_gens_before_sampling
        n_expected_matches_per_species_pair = pop_sizes[1] * n_matches_per_ind[1] + n_expected_samples
        n_expected_total_matches = n_expected_matches_per_species_pair * n_species_pairs
    else
        n_expected_matches_per_species_pair = length(pop_sizes) == 2 ? reduce(*,pop_sizes) : pop_sizes[1]^2
        n_expected_total_matches = sum(pop_sizes[i]*pop_sizes[j]
                                       for i in 1:length(pop_sizes)
                                       for j in (i+1):length(pop_sizes))
    end
    # Create species based on pop_sizes
    g = BasicVectorGenotype([0.0])
    pops, id = [], 1
    for pop_size in pop_sizes
        pop = [BasicIndividual(i, g, Int[]) for i in id:id+pop_size-1]
        push!(pops, pop)
        id += pop_size
    end
    mids = [Int(ps / 2) for ps in pop_sizes]
    species = [PhylogeneticSpecies("species$i", pop[1:mids[i]], pop[mids[i]+1:end]) for (i, pop) in enumerate(pops)]
    rng = StableRNG(1)
    # Create all pairwise matches
    all_matches = []
    for idx1 in 1:length(pop_sizes)-1
        for idx2 in (idx1+1):length(pop_sizes)

            rcmm = RandomCohortMatchMaker(n_matches_per_ind=Dict(
                                    "species$i"=>n_matches_per_ind[i] 
                                    for (i, pop) in enumerate(pops)),
                                  n_samples=n_samples,
                                  gen=gen,
                                  n_gens_before_sampling=n_gens_before_sampling)
            indicies = [idx1, idx2]
            matches = make_matches(rcmm, rng, "interaction1", species[idx1], species[idx2])
            append!(all_matches, matches)
            @test length(matches) == n_expected_matches_per_species_pair
            wrong_num_matches = false
            for (i, (pop, n_matches_for_ind)) in enumerate(zip(pops[indicies], n_matches_per_ind[indicies]))
                wrong_num_matches = false
                for ind in pop
                    num_made_matches_for_ind = length([1 for m in matches if ind.id in m.individual_ids])
                    if num_made_matches_for_ind < n_matches_for_ind
                        wrong_num_matches = true
                        @assert false "individual $(ind.id) made $(num_made_matches_for_ind) matches, expected $(n_matches_for_ind)"
                    end
                end
                @test !wrong_num_matches
            end
        end
    end
    @test length(all_matches) == n_expected_total_matches
    all_matches
end
@testset "RandomCohortMatchMaker" begin
    @testset "10v10, 2 matches per ind" begin
        test_random_cohort_matches([10, 10], [2, 2], gen=2)
        test_random_cohort_matches([20, 20], [4, 4], gen=2)
        # TODO add more tests
    end
    @testset "100v10" begin
        # test all vs all for first gen
        test_random_cohort_matches([100, 10], [1, 10], gen=0, n_samples=10, n_expected_samples=0)
        # test no samples for second gen
        test_random_cohort_matches([100, 10], [1, 10], gen=1, n_gens_before_sampling=10, n_samples=10, n_expected_samples=0)
        test_random_cohort_matches([100, 10], [1, 10], gen=2, n_samples=10, n_expected_samples=10)
        test_random_cohort_matches([100, 10], [2, 20], gen=2)
        test_random_cohort_matches([100, 10], [5, 50], gen=2)
        test_random_cohort_matches([100, 10], [10, 100], gen=2)
    end
    @testset "3SpeciesSameSizeSameNumMatches" begin
        test_random_cohort_matches([100,100,100], [1,1,1], gen=0, n_samples=10, n_expected_samples=0)
        test_random_cohort_matches([100,100,100], [1,1,1], gen=2, n_samples=10, n_expected_samples=10)
        test_random_cohort_matches([100,100,100], [10,10,10], gen=2, n_samples=10, n_expected_samples=10)
        test_random_cohort_matches([100,100,100], [100,100,100], gen=2)
    end
end

function test_parent_v_children_matches(n_parents::Vector{Int},
                                        n_children::Vector{Int};
                                        n_samples::Int,
                                        all_vs_all_expected::Bool=false
                                        )
    "Function goes over each species pair and tests that the correct number of matches are made.
    The number of matches made between two species is the same for all species pairs.
    The number of matches made by an individual is the same for all individuals in a species.
    Arguments:
        n_samples::Int: Number of additional child v child estimates to compute
    "
    @assert length(n_parents) == length(n_children)
    n_species_pairs = length(n_parents) * (length(n_parents) - 1) / 2 

    if all_vs_all_expected
        n_expected_total_matches = 0
        for idx1 in 1:length(n_parents)-1
            for idx2 in (idx1+1):length(n_parents)
                n_expected_total_matches += (n_parents[idx1] + n_children[idx1]) * (n_parents[idx2] + n_children[idx2])
            end
        end
    else
        # pairwise matches between species
        n_expected_total_matches = 0
        # Count parents vs parents
        for i in 1:length(n_parents)-1, j in (i+1):length(n_parents)
            n_expected_total_matches += n_parents[i] * n_parents[j]
        end
        # Count parents vs children
        for i in 1:length(n_parents), j in [1:(i-1);(i+1):length(n_parents)]
            n_expected_total_matches += n_parents[i] * n_children[j]
        end
        n_expected_total_matches += n_samples * n_species_pairs
    end
    # Create species based on n_parents and n_children
    species = make_dummy_phylo_species(n_parents, n_children, first_gen=all_vs_all_expected)
    rng = StableRNG(1)
    pvcmm = ParentsVsAllMatchMaker(n_samples=n_samples)

    # Create all parent v child matches across all species
    all_matches = []
    for idx1 in 1:length(species)-1
        for idx2 in (idx1+1):length(species)
            matches = make_matches(pvcmm, rng, "interaction1", species[idx1], species[idx2])
            append!(all_matches, matches)

            if all_vs_all_expected
                @test length(matches) == (n_parents[idx1] + n_children[idx1]) * (n_parents[idx2] + n_children[idx2])
                continue
            end
            @test length(matches) ==
                  (n_parents[idx1] * n_parents[idx2]) +
                  (n_parents[idx1] * n_children[idx2]) +
                  (n_parents[idx2] * n_children[idx1]) +
                  n_samples
            wrong_num_matches = false
            n_parent_matches = 0
            for (i, s) in enumerate(species[[idx1,idx2]])
                wrong_num_matches = false
                # assert each parent is matched up with the entire pop
                for ind in s.population
                    num_made_matches_for_p = length([1 for m in matches if ind.id in m.individual_ids])
                    n_parent_matches += num_made_matches_for_p
                    if num_made_matches_for_p != length(s.children) + length(s.population)
                        wrong_num_matches = true
                        @assert false "individual $(ind.id) made $(num_made_matches_for_p) matches, expected $(length(s.children))"
                    end
                end
                @test !wrong_num_matches
            end
            # This code counts parents vs parents twice, so we add them back in
            @test length(matches) + (n_parents[idx1] * n_parents[idx2]) == n_parent_matches + n_samples
        end
    end
    @test length(all_matches) == n_expected_total_matches
    all_matches
end


@testset "ParentsVsAllMatchMaker" begin
    test_parent_v_children_matches([10, 10], [10, 10], n_samples=1)
    test_parent_v_children_matches([10, 10], [10, 10], n_samples=1, all_vs_all_expected=true)
    test_parent_v_children_matches([10, 10], [100, 100], n_samples=100)
    test_parent_v_children_matches([10, 10, 10], [10, 10, 10], n_samples=1)
    test_parent_v_children_matches([10, 10, 10], [10, 10, 10], n_samples=1, all_vs_all_expected=true)
end
