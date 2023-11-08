using CoEvo.Individuals.Basic: BasicIndividual
using CoEvo.Genotypes.Vectors: BasicVectorGenotype
function test_random_cohort_matches(pop_sizes::Vector{Int},
                                    n_matches_per_ind::Vector{Int};
                                    include_children::Bool=false)
    "Function goes over each species pair and tests that the correct number of matches are made.
    The number of matches made between two species is the same for all species pairs.
    The number of matches made by an individual is the same for all individuals in a species.
    Arguments:
        include_children::Bool: If true, split pop between pop and children, otherwise only use pop 
    "
    # With random cohort matches, the number of matches made between two species is the same
    # for all species pairs.
    @assert length(Set(ps * nm for (ps, nm) in zip(pop_sizes, n_matches_per_ind))) == 1
    n_expected_matches_per_species_pair = pop_sizes[1] * n_matches_per_ind[1]
    # sum(1:length(pop_sizes)-1) is the number of species pairs, each pair has the same num of matches
    n_expected_total_matches =n_expected_matches_per_species_pair * sum(1:length(pop_sizes)-1)
    # Create species based on pop_sizes
    g = BasicVectorGenotype([0.0])
    pops, id = [], 1
    for pop_size in pop_sizes
        pop = [BasicIndividual(i, g, Int[]) for i in id:id+pop_size-1]
        push!(pops, pop)
        id += pop_size
    end
    if include_children
        mids = [Int(ps / 2) for ps in pop_sizes]
        species = [PhylogeneticSpecies("species$i", pop[1:mids[i]], pop[mids[i]+1:end]) for (i, pop) in enumerate(pops)]
    else
        species = [PhylogeneticSpecies("species$i", pop, typeof(pop)()) for (i, pop) in enumerate(pops)]
    end
    rng = StableRNG(1)
    rcmm = RandomCohortMatchMaker(n_matches_per_ind=Dict("species$i"=>n_matches_per_ind[i] for (i, pop) in enumerate(pops))) 
    # Create all pairwise matches
    all_matches = []
    for idx1 in 1:length(pop_sizes)-1
        for idx2 in (idx1+1):length(pop_sizes)
            indicies = [idx1, idx2]
            matches = make_matches(rcmm, rng, "interaction1", species[idx1], species[idx2])
            append!(all_matches, matches)
            @test length(matches) == n_expected_matches_per_species_pair
            wrong_num_matches = false
            for (i, (pop, n_matches_for_ind)) in enumerate(zip(pops[indicies], n_matches_per_ind[indicies]))
                wrong_num_matches = false
                for ind in pop
                    num_made_matches_for_ind = length([1 for m in matches if ind.id in m.individual_ids])
                    if num_made_matches_for_ind != n_matches_for_ind
                        wrong_num_matches = true
                        @assert num_made_matches_for_ind == n_matches_for_ind "individual $(ind.id) made $(num_made_matches) matches, expected $(n_matches_for_ind)"
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
        test_random_cohort_matches([10, 10], [2, 2])
        test_random_cohort_matches([20, 20], [4, 4], include_children=true)
    end
    @testset "100v10" begin
        test_random_cohort_matches([100, 10], [1, 10])
        test_random_cohort_matches([100, 10], [2, 20], include_children=true)
        test_random_cohort_matches([100, 10], [5, 50])
        test_random_cohort_matches([100, 10], [10, 100], include_children=true)
    end
    @testset "3SpeciesSameSizeSameNumMatches" begin
        test_random_cohort_matches([100,100,100], [1,1,1])
        test_random_cohort_matches([100,100,100], [10,10,10], include_children=true)
        test_random_cohort_matches([100,100,100], [100,100,100])
    end
end
