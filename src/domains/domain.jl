module Domains
include("./numbers_game/numbers_game.jl")
using .NumbersGame
include("./sorting-network/domain.jl")
using .SortingNetwork
end
