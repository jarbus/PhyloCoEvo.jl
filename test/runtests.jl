using PhyloCoEvo
using StableRNGs
using CoEvo
using Test

@testset "PhyloCoEvo.jl" begin
    # Write your tests here.
    include("x/ng/1/config.jl")
    include("sorting-networks.jl")
    include("x/sn/1/config.jl")
    include("x/pg/1/config.jl")
    include("estimation.jl")
    include("matchmakers.jl")
    include("x/ng-mm/1/config.jl")
end
