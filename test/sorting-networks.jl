using CoEvo.Counters.Basic: BasicCounter
using CoEvo.Environments.Stateless: StatelessEnvironment
using PhyloCoEvo.Genotypes: create_genotypes
using PhyloCoEvo.Phenotypes: create_phenotype
using PhyloCoEvo: SortingNetworkDomain
@testset "SortingNetworks" begin
    # Test that everything runs
    sngc = SortingNetworkGenotypeCreator(1, 2)
    snpc = SortingNetworkPhenotypeCreator(2)
    # create sorting network genotypes and phenotypes
    rng = StableRNG(1)
    gene_id_counter = BasicCounter(0)
    n_pop = 1
    sngenotypes = create_genotypes(sngc, rng, gene_id_counter, n_pop)
    @test length(sngenotypes) == 1
    phenotype = create_phenotype(snpc, sngenotypes[1])
    # create sorting network domain
    domain = SortingNetworkDomain(Partial())


    @testset "Genotype" begin
        # no activate comparators
        Codon = PhyloCoEvo.SortingNetworkCodon
        codons = (Codon(1, 0b0000000000100001),)
        genotype = SortingNetworkGenotype(codons, 2)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 2
        @test phenotype.network == zeros(Int64, 0, 2)
                                    
        # one active comparator, encoding 2 1
        codons = (Codon(1, 0b1111111100010000),)
        genotype = SortingNetworkGenotype(codons, 2)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 2
        @test phenotype.network == [2 1]

        # two active comparators encoding 15 7, 3 1,
        # one inactive comparator encoding 7 3
        codons = (Codon(1, 0b1111101111110111),
                  Codon(2, 0b1101111100110001),
                  Codon(3, 0b0000001101110011))
        genotype = SortingNetworkGenotype(codons, 16)
        snpc = SortingNetworkPhenotypeCreator(16)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 16
        @test phenotype.network == [16 8; 4 2]
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

            codons = PhyloCoEvo.SortingNetworkCodon[]
            c = 1
            for i in 1:n-1
                for j in i:-1:1
                    codon = 0b1111111100000000
                    codon = codon | ((j-1) << 4)
                    codon = codon | (j << 0)
                    codon = PhyloCoEvo.SortingNetworkCodon(c, codon)
                    push!(codons, codon)
                end
            end
            SortingNetworkGenotype(Tuple(codons), n)
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
            SortingNetworkPhenotype(network, n)
        end

        @testset "CorrectSorting" begin
            # Sanity check
            snp = SortingNetworkPhenotype([1 2;], 2)
            @test [1,2] == PhyloCoEvo.netsort(snp, (1, 2))
            @test [1,2] == PhyloCoEvo.netsort(snp, (2, 1))

            # Wikipedia example, 4 inputs
            # __________
            # _|___|____
            # _|_|___|__
            # ___|_|____
            snp = SortingNetworkPhenotype([1 3; 2 4; 1 2; 3 4; 2 3;], 4)
            sorted = [1, 2, 3, 4]
            for perm in permute(sorted)
                @assert sorted == PhyloCoEvo.netsort(snp, Tuple(perm))
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
                @assert sorted == PhyloCoEvo.netsort(snp, Tuple(perm))
                @assert sorted == PhyloCoEvo.netsort(sngp, Tuple(perm))
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
                    @assert sorted == PhyloCoEvo.netsort(snp, Tuple(bits))
                    @assert sorted == PhyloCoEvo.netsort(sngp, Tuple(bits))
                end
            end
            test_size(16)
            @test true

            # Test percent_sorted
            @test PhyloCoEvo.percent_sorted(sng) == 1.0
            sng_missing = SortingNetworkGenotype(sng.codons[1:5], 16)
            @test PhyloCoEvo.percent_sorted(sng_missing) < 1.0
        end

        # TODO: test that wires are valid with smaller networks

        @testset "IncorrectSorting" begin
            # Empty network
            snp = SortingNetworkPhenotype(zeros(Int64,0,2), 2)
            @test [2,1] == PhyloCoEvo.netsort(snp, (2, 1))

            # Four inputs, one comparator
            snp = SortingNetworkPhenotype([1 2;], 4)
            @test [3,4,2,1] == PhyloCoEvo.netsort(snp, (4,3,2,1))
        end
    end


    @testset "Mutator" begin
        snm = SortingNetworkMutator()
        # test one inactive comparator
        Codon = PhyloCoEvo.SortingNetworkCodon
        codons = (Codon(1, 0b0000000000100001),)
        genotype = SortingNetworkGenotype(codons, 1)
        new_geno = mutate(snm, rng, gene_id_counter, genotype)
        # test that one of the first 8 activation bits is flipped
        @test count_ones(new_geno.codons[1].data & 0b1111111100000000) == 1
        # test that one of the last 8 bit flip bits is flipped
        @test new_geno.codons[1].data & 0b0000000011111111 == 0b0000000000100001

        # test one active comparator
        codons = (Codon(1, 0b1111111100100001),)
        genotype = SortingNetworkGenotype(codons, 1)
        new_geno = mutate(snm, rng, gene_id_counter, genotype)
        # test that one of the first 8 activation bits is flipped (8 zeros at the end)
        @test count_zeros(new_geno.codons[1].data & 0b1111111100000000) == 9
        # test that one of the last 8 bit flip bits is flipped
        @test count_ones(new_geno.codons[1].data & 0b0000000011111111) ∈ [1, 3]
    end
end

@testset "SortingNetworkTestCase" begin
    # create sorting network test case genotypes and phenotypes
    rng = StableRNG(1)
    gene_id_counter = BasicCounter(0)
    n_pop = 1
    sntcgc = SortingNetworkTestCaseGenotypeCreator(16)
    sntcpc = SortingNetworkTestCasePhenotypeCreator(16)
    sntcgenotypes = create_genotypes(sntcgc, rng, gene_id_counter, n_pop)
    @test length(sntcgenotypes) == 1
    @test sum([sntcgenotypes[1].inputs...] .== 1:16) ∈ 12:16
    @test sort([sntcgenotypes[1].inputs...]) == 1:16
    phenotype = create_phenotype(sntcpc, sntcgenotypes[1])

    x = PhyloCoEvo.swap(rng,4,0)
    @test x == 1:4

    x = PhyloCoEvo.swap(rng,4,1)
    @test length(x) == 4
    @test sum(x .== 1:4) == 2

    x = PhyloCoEvo.swap(rng,16,2)
    @test length(x) == 16
    @test sum(x .== 1:16) ∈ 12:16

    x = PhyloCoEvo.swap(rng,16,100)
    @test length(x) == 16
    @test sum(x .== 1:16) ∉ 14:16

    # Test that initial test case is nearly sorted
    sntcgenotypes = create_genotypes(sntcgc, rng, gene_id_counter, 5)
    for sntcg in sntcgenotypes
        @test sum([sntcg.inputs...] .== 1:16) ∈ 12:16
    end


    @testset "Mutator" begin
        sntcm = SortingNetworkTestCaseMutator()
        # test one swap
        genotype = SortingNetworkTestCaseGenotype(1, (1,2,3,4))
        new_geno = mutate(sntcm, rng, gene_id_counter, genotype)
        # test that two inputs are different than 1:4
        @test sum([new_geno.inputs...] .== 1:4) == 2
        @test sort([new_geno.inputs...]) == 1:4
        # test two swaps
        genotype = SortingNetworkTestCaseGenotype(1, (1,2,3,4))
        new_geno = mutate(sntcm, rng, gene_id_counter, genotype)
        new_geno = mutate(sntcm, rng, gene_id_counter, new_geno)
        @test sort([new_geno.inputs...]) == 1:4
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
    snp = SortingNetworkPhenotype([1 3; 2 4; 1 2; 3 4; 2 3;], 4)
    sntc = SortingNetworkTestCasePhenotype((4,3,2,1))
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [4.0, 0.0]

    # Test that a network with no comparators results in 
    # zero fitness in the outcome when tested against an incorrect list
    snp = SortingNetworkPhenotype(zeros(Int,0,2), 4)
    sntc = SortingNetworkTestCasePhenotype((4,3,2,1))
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [0.0, 4.0]


    # Test that a network which sorts two numbers correctly
    # and two numbers incorrectly results in a fitness of 2.0
    # for network and inputs
    snp = SortingNetworkPhenotype([1 4; ], 4)
    sntc = SortingNetworkTestCasePhenotype((4,3,2,1))
    env = StatelessEnvironment(domain, [snp, sntc])
    outcome = PhyloCoEvo.Environments.get_outcome_set(env)
    @test outcome == [2.0, 2.0]
end
