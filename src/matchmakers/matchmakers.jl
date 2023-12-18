module MatchMakers

include("./random_cohort/random_cohort.jl")
using .RandomCohort

include("./parents_vs_all/parents_vs_all.jl")
using .ParentsVsAll

end

