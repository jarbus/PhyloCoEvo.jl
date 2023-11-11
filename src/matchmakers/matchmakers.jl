module MatchMakers

include("./random_cohort/random_cohort.jl")
using .RandomCohort

include("./parents_vs_children/parents_vs_children.jl")
using .ParentsVsChildren

end

