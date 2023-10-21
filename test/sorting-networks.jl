@testset "SortingNetworks" begin
    sngc = SortingNetworkGenotypeCreator(1, 1)
    snpc = SortingNetworkPhenotypeCreator(1)
    rng = StableRNG(1)
    gene_id_counter = Counter(0)
    n_pop = 1
    genotypes = create_genotypes(sngc, rng, gene_id_counter, n_pop)
    @test length(genotypes) == 1
    phenotype = create_phenotype(snpc, genotypes[1])
    domain = SortingNetworkDomain(Partial())

    @testset "Genotype" begin
        # no activate comparators
        Codon = PhyloCoEvo.SortingNetworkCodon
        codons = (Codon(1, 0b0000000000100001),)
        genotype = SortingNetworkGenotype(codons, 1)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 1
        @test phenotype.network == zeros(Int64, 0, 2)
                                    
        # one active comparator, encoding 2 1
        codons = (Codon(1, 0b1111111100100001),)
        genotype = SortingNetworkGenotype(codons, 1)
        phenotype = create_phenotype(snpc, genotype)
        @test phenotype.n == 1
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
        @test phenotype.network == [15 7; 3 1]
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

        function create_insertion_sort_network(n::Int)
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
            snp = create_insertion_sort_network(4)
            sorted = [1, 2, 3, 4]
            for perm in permute(sorted)
                @assert sorted == PhyloCoEvo.netsort(snp, Tuple(perm))
            end
            @test true

            # Test all possible combinations of zeros and ones are sorted correctly for size 16
            function test_size(n)
                snp = create_insertion_sort_network(n)
                for i in 0:2^n
                    # create a 16-bit vector of zeros and ones of i
                    bits = [Int(d) for d in digits(i, base=2, pad=n)]
                    sorted = sort(bits)
                    @assert sorted == PhyloCoEvo.netsort(snp, Tuple(bits))
                end
            end
            test_size(16)
            @test true
        end

        @testset "IncorrectSorting" begin
            # Empty network
            snp = SortingNetworkPhenotype(zeros(Int64,0,2), 2)
            @test [2,1] == PhyloCoEvo.netsort(snp, (2, 1))

            # Four inputs, one comparator
            snp = SortingNetworkPhenotype([1 2;], 4)
            @test [3,4,2,1] == PhyloCoEvo.netsort(snp, (4,3,2,1))
        end
    end
end
