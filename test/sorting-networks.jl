using CoEvo.Counters.Basic: BasicCounter
using CoEvo.Environments.Stateless: StatelessEnvironment
using PhyloCoEvo.Genotypes.SortingNetwork: SortingNetworkGenotype, SortingNetworkTestCaseGenotype, SortingNetworkCodon, SortingNetworkGenotypeCreator, SortingNetworkTestCaseGenotypeCreator, random_codon, swap
using PhyloCoEvo.Phenotypes.SortingNetwork: SortingNetworkPhenotype, SortingNetworkTestCasePhenotype, netsort, SortingNetworkPhenotypeCreator, SortingNetworkTestCasePhenotypeCreator
using PhyloCoEvo.Metrics.SortingNetwork: percent_sorted
using PhyloCoEvo.Mutators.SortingNetwork: SortingNetworkMutator, SortingNetworkTestCaseMutator
using PhyloCoEvo.Domains.SortingNetwork: SortingNetworkDomain, Partial
using CoEvo.Phenotypes: create_phenotype
using CoEvo.Genotypes: create_genotypes

@testset "SortingNetworks" begin
    # Test that everything runs
    sngc = SortingNetworkGenotypeCreator(1, 2, 0, 100)
    snpc = SortingNetworkPhenotypeCreator(2)
    # create sorting network genotypes and phenotypes
    rng = StableRNG(1)
    gene_id_counter = BasicCounter(0)
    n_pop = 1
    sngenotypes = create_genotypes(sngc, rng, gene_id_counter, n_pop)
    @test length(sngenotypes) == 1
    phenotype = create_phenotype(snpc, sngenotypes[1])
    snm = SortingNetworkMutator()
    new_geno = mutate(snm, rng, gene_id_counter, sngenotypes[1])
    # create sorting network domain
    domain = SortingNetworkDomain(Partial())


    @testset "Genotype" begin
        # no activate comparators
        Codon = SortingNetworkCodon
        codons = Codon[]
        genotype = SortingNetworkGenotype(codons, 2, 0, 100)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 2
        @test phenotype.network == zeros(Int64, 0, 2)
     
        # one active comparator, encoding 2 1
        codons = [Codon(1, 2, 1)]
        genotype = SortingNetworkGenotype(codons, 2, 0, 100)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 2
        @test phenotype.network == [2 1]

        # Create phenotype from genotype
        codons = [Codon(1, 16, 8),
                  Codon(2, 4, 2)]
        genotype = SortingNetworkGenotype(codons, 16, 0, 100)
        snpc = SortingNetworkPhenotypeCreator(16)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 16
        @test phenotype.network == [16 8; 4 2]

        # Create genotypes 80 codons 16 inputs
        # Test that swaps are seeded
        n_pop16 = 10
        sngc = SortingNetworkGenotypeCreator(80, 16, 60, 100)
        sngenotypes = create_genotypes(sngc, rng, gene_id_counter, n_pop)
        @test length(sngenotypes[1].codons) == 80
        snpc16 = SortingNetworkPhenotypeCreator(16)
        swaps = [ 
            # Column 1: swap adjacent
            1  2; 3  4; 5  6; 7  8; 9 10; 11 12; 13 14; 15 16;
            # Columns 2, 3: swap 1 apart
            1 3; 5 7; 9 11; 13 15;
            2 4; 6 8; 10 12; 14 16;
            # Columns 4-7: Swap 3 apart
            1 5; 9 13;
            2 6; 10 14;
            3 7; 11 15;
            4 8; 12 16;
            # Columns 8-15: Swap 7 apart
            1 9; 2 10; 3 11; 4 12; 5 13; 6 14; 7 15; 8 16;
        ]
        n_seed_swaps = size(swaps, 1)
        for sng in sngenotypes
            snphenotype = create_phenotype(snpc16, sng)
            @test snphenotype.network[1:n_seed_swaps, :] == swaps
        end
    end
    
    @testset "Phenotype" begin
        function permute(arr)
            length(arr) == 0 && return [Int[]]
            perms = []
            for i in 1:length(arr)
                rest = permute(vcat(arr[1:i-1], arr[i+1:end]))
                for p in rest
                    push!(perms, vcat(arr[i], p))
                end
            end
            return perms
        end
    
        function create_insertion_sort_genotype(n::Int)
            num_codons = Int(n * (n - 1) / 2)
            codons = SortingNetworkCodon[]
            c = 1
            for i in 1:n-1
                for j in i:-1:1
                    codon = SortingNetworkCodon(c, j, j+1)
                    push!(codons, codon)
                end
            end
            SortingNetworkGenotype(codons, n, num_codons-1, num_codons+1)
        end
    
        function create_insertion_sort_phenotype(n::Int)
            num_comparators = Int(n * (n - 1) / 2)
            network = zeros(Int64, num_comparators, 2)
            comparator = 1
            for i in 1:n-1
                for j in i:-1:1
                    network[comparator, :] = [j, j+1]
                    comparator += 1
                end
            end
            SortingNetworkPhenotype(network, n,n-1, n+1)
        end
    
        @testset "CorrectSorting" begin
            # Sanity check
            snp = SortingNetworkPhenotype([1 2;], 2, 0, 2)
            @test [1,2] == netsort(snp, [1, 2])
            @test [1,2] == netsort(snp, [2, 1])
    
            # Wikipedia example, 4 inputs
            # __________
            # _|___|____
            # _|_|___|__
            # ___|_|____
            snp = SortingNetworkPhenotype([1 3; 2 4; 1 2; 3 4; 2 3;], 4, 0, 10)
            sorted = [1, 2, 3, 4]
            for perm in permute(sorted)
                @assert sorted == netsort(snp, perm)
            end
            @test true
    
            # Insertion sort example, n inputs
            # Example: for n = 4
            # __________
            # _|_|_|___
            # __|_|___
            # ___|_____
            sng = create_insertion_sort_genotype(4) # test raw phenotype
            sngp = create_phenotype(snpc, sng)      # test genotype->phenotype
            snp = create_insertion_sort_phenotype(4)
            sorted = [1, 2, 3, 4]
            for perm in permute(sorted)
                @assert sorted == netsort(snp, perm)
                @assert sorted == netsort(sngp, perm)
            end
            @test true
    
            # Test all possible combinations of zeros and ones are sorted correctly for size 16
            function test_size(n)
                snp = create_insertion_sort_phenotype(n)
                snpc = SortingNetworkPhenotypeCreator(n)
                sng = create_insertion_sort_genotype(n)
                sngp = create_phenotype(snpc, sng)
                for i in 0:2^n
                    # create a 16-bit vector of zeros and ones of i
                    bits = [Int(d) for d in digits(i, base=2, pad=n)]
                    sorted = sort(bits)
                    @assert sorted == netsort(snp, bits)
                    @assert sorted == netsort(sngp, bits)
                end
            end
            test_size(16)
            @test true
    
            # Test percent_sorted
            @test percent_sorted(sng) == 1.0
            sng_missing = SortingNetworkGenotype(sng.codons[1:5], 16, 0, 100)
            @test percent_sorted(sng_missing) < 1.0
        end
    
        # TODO: test that wires are valid with smaller networks
    
        @testset "IncorrectSorting" begin
            # Empty network
            snp = SortingNetworkPhenotype(zeros(Int64,0,2), 2, 0, 2)
            @test [2,1] == netsort(snp, [2, 1])
    
            # Four inputs, one comparator
            snp = SortingNetworkPhenotype([1 2;], 4, 0, 4)
            @test [3,4,2,1] == netsort(snp, [4,3,2,1])
        end
    end


    @testset "Mutator" begin
        Codon = SortingNetworkCodon
        snm = SortingNetworkMutator()
        # test one active comparator
        codons = [Codon(1, 1, 2),]
        genotype = SortingNetworkGenotype(codons, 16, 1, 2)

        # Test move
        mrng = StableRNG(1)
        new_geno = mutate(snm, mrng, gene_id_counter, genotype)
        @test length(new_geno.codons) == 1

        # Test delete
        codons = [Codon(1, 1, 2), Codon(1, 1, 2)]
        genotype = SortingNetworkGenotype(codons, 16, 1, 2)
        mrng = StableRNG(2)
        new_geno = mutate(snm, mrng, gene_id_counter, genotype)
        @test length(new_geno.codons) == 1

        # Test insert
        codons = [Codon(1, 1, 2),]
        genotype = SortingNetworkGenotype(codons, 16, 1,2)
        mrng = StableRNG(3)
        new_geno = mutate(snm, mrng, gene_id_counter, genotype)
        @test length(new_geno.codons) == 2

        # Test rewire
        mrng = StableRNG(8)
        codons = [Codon(1, 1, 2),]
        genotype = SortingNetworkGenotype(codons, 16, 1,2)
        new_geno = mutate(snm, mrng, gene_id_counter, genotype)
        @test length(new_geno.codons) == 1
        @test new_geno.codons[1].one != 1 || new_geno.codons[1].two != 2
    end
end

@testset "SortingNetworkTestCase" begin
    # create sorting network test case genotypes and phenotypes
    rng = StableRNG(1)
    gene_id_counter = BasicCounter(0)
    n_pop = 1
    sntcgc = SortingNetworkTestCaseGenotypeCreator(20, 16)
    sntcpc = SortingNetworkTestCasePhenotypeCreator(16)
    sntcgenotypes = create_genotypes(sntcgc, rng, gene_id_counter, n_pop)
    @test length(sntcgenotypes) == 1
    @test sum([sntcgenotypes[1].tests[1]...] .== 1:16) ∈ 12:16
    @test sort([sntcgenotypes[1].tests[1]...]) == 1:16
    phenotype = create_phenotype(sntcpc, sntcgenotypes[1])

    x = swap(rng,4,0)
    @test x == 1:4

    x = swap(rng,4,1)
    @test length(x) == 4
    @test sum(x .== 1:4) == 2

    x = swap(rng,16,2)
    @test length(x) == 16
    @test sum(x .== 1:16) ∈ 12:16

    x = swap(rng,16,100)
    @test length(x) == 16
    @test sum(x .== 1:16) ∉ 14:16

    # Test that initial test case is nearly sorted
    sntcgenotypes = create_genotypes(sntcgc, rng, gene_id_counter, 5)
    for sntcg in sntcgenotypes
        @test sum([sntcg.tests[1]...] .== 1:16) ∈ 12:16
    end


    @testset "Mutator" begin
        sntcm = SortingNetworkTestCaseMutator()
        # test one swap
        genotype = SortingNetworkTestCaseGenotype(1, [[1,2,3,4]])
        new_geno = mutate(sntcm, rng, gene_id_counter, genotype)
        # test that two inputs are different than 1:4
        @test sum([new_geno.tests[1]...] .== 1:4) == 2
        @test sort([new_geno.tests[1]...]) == 1:4
        # test two swaps
        genotype = SortingNetworkTestCaseGenotype(1, [[1,2,3,4]])
        new_geno = mutate(sntcm, rng, gene_id_counter, genotype)
        new_geno = mutate(sntcm, rng, gene_id_counter, new_geno)
        @test sort([new_geno.tests[1]...]) == 1:4
    end
end


@testset "SortingNetworkDomain" begin

    domain = SortingNetworkDomain(Partial())

    # Wikipedia example, correct, 4 inputs
    # __________
    # _|___|____
    # _|_|___|__
    # ___|_|____
    # Test that a network with comparators results in 
    # full fitness in the outcome when tested against a correct list
    snp = SortingNetworkPhenotype([1 3; 2 4; 1 2; 3 4; 2 3;], 4, 5, 5)
    sntc = SortingNetworkTestCasePhenotype([[4,3,2,1], [1,2,3,4], [2,3,1,4]])
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [3.0, 0.0]

    # Test that a correct maximum-length network produces a outcome of 3
    snp_extra = SortingNetworkPhenotype([1 3; 2 4; 1 2; 3 4; 2 3; 1 4;], 4, 5, 6)
    sntc = SortingNetworkTestCasePhenotype([[4,3,2,1], [1,2,3,4], [2,3,1,4]])
    env = StatelessEnvironment(domain, [snp_extra, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [3.0, 0.0]

    # Test that a network with no comparators results in 
    # zero fitness in the outcome when tested against an incorrect list
    snp = SortingNetworkPhenotype(zeros(Int,0,2), 4, 0, 0)
    sntc = SortingNetworkTestCasePhenotype([[4,3,2,1], [2,3,1,4], [3,4,1,2]])
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [0.0, 3.0]


    # Test that a network which sorts the first and last numbers correctly
    # can get some cases right and some cases wrong
    snp = SortingNetworkPhenotype([1 4; ], 4, 1, 1)
    sntc = SortingNetworkTestCasePhenotype([[4,3,2,1], [4,2,3,1], [2,3,1,4], [1,2,3,4]])
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [2.0, 2.0]

    # Test that a network which gets no test cases right gets 0
    snp = SortingNetworkPhenotype([1 4; ], 4, 0, 1)
    sntc = SortingNetworkTestCasePhenotype([[4,3,2,1], [4,3,2,1], [4,3,2,1]])
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [0.0, 3.0]

end
