using PhyloCoEvo
using StableRNGs
using CoEvo
using Test

@testset "PhyloCoEvo.jl" begin
    # Write your tests here.
    #include("x/ng/1/config.jl")
    include("sorting-networks.jl")
    include("x/sn/config.jl")
end
