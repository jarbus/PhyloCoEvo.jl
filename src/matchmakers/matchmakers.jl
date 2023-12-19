module MatchMakers

include("./random_cohort/random_cohort.jl")
using .RandomCohort

include("./parents_vs_all/parents_vs_all.jl")
using .ParentsVsAll

include("./random_parents_vs_children/random_parents_vs_children.jl")
using .RandomParentsVsChildren

end

