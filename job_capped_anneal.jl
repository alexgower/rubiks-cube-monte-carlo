using Distributed
addprocs(2)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/capped_anneal_experiment.jl")


starting_configuration_L_5 = [[2 3 3 2 6; 5 3 6 4 4; 1 2 5 1 5; 2 2 6 6 2; 2 1 2 3 2],[1 2 5 5 3; 6 6 1 2 3; 2 5 5 4 3; 4 1 6 5 4; 4 6 6 2 4], [3 4 5 4 6; 2 5 4 4 1; 5 3 2 5 4; 4 6 1 2 5; 2 5 5 4 6], [3 3 1 6 3; 2 3 1 6 3; 6 6 4 3 2; 1 1 2 3 2; 6 3 1 5 3], [2 3 4 4 1; 3 3 1 3 6; 3 4 4 5 1; 4 5 5 4 1; 6 1 1 6 6], [4 1 2 3 5; 2 1 6 6 5; 4 1 5 5 1; 2 1 5 1 6; 5 4 4 3 6]]
glassy_configuration_L_5 = [[5 3 4 4 3; 3 3 4 6 4; 5 5 5 6 6; 1 5 5 6 6; 1 2 4 4 4], [5 5 5 2 6; 5 5 2 2 1; 5 5 5 1 1; 5 5 1 1 6; 6 2 1 3 6], [3 3 1 1 6; 5 3 1 1 1; 2 5 2 1 1; 2 4 4 4 4; 2 4 4 4 4], [2 2 2 2 1; 2 2 2 2 2; 3 6 4 4 2; 6 6 6 6 3; 6 6 6 3 3], [3 6 5 5 3; 1 1 6 3 3; 1 1 4 3 3; 4 1 6 6 6; 6 2 2 2 2], [2 3 5 4 4; 3 3 5 4 4; 3 3 5 1 1; 3 3 5 1 1; 2 5 4 4 6]] 


starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 4 1 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 6 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 2; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 6 5 2 2; 3 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 4 3 1 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 2 1 6 3; 6 3 6 6 5 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
glassy_configuration_L_7 = [[6 4 4 2 2 2 2; 4 4 4 4 4 4 4; 3 2 6 4 4 4 4; 3 3 5 4 4 4 4; 2 2 5 4 4 4 1; 5 5 5 1 5 4 3; 1 5 6 6 5 5 5], [1 1 3 3 5 4 4; 6 6 6 6 6 5 4; 6 6 6 6 6 6 4; 6 6 6 6 1 3 3; 6 6 4 4 5 3 3; 1 6 1 3 3 3 3; 2 2 3 3 3 6 6], [6 6 4 4 6 5 3; 4 4 4 4 4 5 4; 3 3 6 6 6 4 4; 3 3 3 6 3 1 5; 3 3 3 3 2 1 5; 6 5 3 3 3 5 1; 2 3 4 4 4 6 5], [5 5 5 2 3 5 1; 5 5 5 5 1 1 1; 5 5 6 1 1 1 1; 1 1 1 1 5 2 1; 5 5 5 5 5 5 3; 2 5 2 5 5 5 5; 2 2 2 2 1 5 5], [2 2 2 2 4 5 5; 2 2 6 6 6 3 1; 1 2 2 2 6 3 6; 2 2 2 2 1 1 1; 2 2 6 2 1 2 3; 2 6 6 3 3 3 6; 2 3 3 3 3 3 3], [5 2 2 2 4 4 6; 1 1 2 2 4 3 4; 2 6 4 1 4 3 2; 2 1 1 1 1 1 2; 6 3 1 1 1 2 2; 6 1 1 1 1 5 4; 5 1 1 4 6 1 6]]

new_starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 1 4 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 2 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 6; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 3 5 2 2; 6 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 1 3 1 4 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 5 1 6 3; 6 3 6 6 2 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 1 2 4 6 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]


# starting_configuration_L_9 = [[2 4 2 6 4 5 2 6 2; 6 6 3 5 1 3 4 3 4; 5 6 6 3 1 5 1 4 4; 4 1 5 6 6 3 1 4 1; 2 6 4 4 6 2 4 5 5; 5 5 1 2 2 2 1 4 3; 3 4 1 1 5 3 3 3 6; 3 4 2 3 2 4 3 4 2; 6 6 5 1 1 3 5 4 4], [6 6 6 3 2 6 1 2 4; 5 1 5 1 2 6 6 2 4; 4 5 2 3 3 1 5 1 2; 4 4 3 4 6 1 6 5 1; 5 5 4 6 6 3 4 3 6; 5 5 4 6 2 3 4 1 4; 6 6 6 2 4 6 1 6 6; 5 5 1 5 5 1 5 4 1; 4 3 5 5 4 1 4 2 5], [5 2 5 2 6 3 1 2 2; 5 3 3 4 3 5 5 4 1; 3 5 4 2 3 1 6 5 2; 4 5 3 2 5 4 2 2 3; 5 4 5 5 1 3 2 6 1; 6 4 3 6 5 1 3 1 6; 5 4 6 3 1 5 3 2 1; 6 5 1 4 3 1 3 6 3; 2 4 5 6 2 1 1 3 5], [3 2 3 5 4 6 2 6 3; 4 6 3 4 3 2 3 1 5; 4 1 2 1 2 4 2 5 3; 6 6 6 1 1 2 6 6 5; 2 2 6 6 6 2 4 1 3; 1 1 2 1 4 3 2 1 6; 3 2 1 1 4 1 4 3 2; 5 3 2 3 1 3 2 1 3; 5 5 1 6 1 5 5 1 2], [5 3 4 6 4 2 5 5 6; 4 2 5 5 1 1 5 2 1; 6 1 2 4 6 2 1 1 4; 6 3 1 5 1 6 2 5 6; 1 2 2 3 3 2 5 5 4; 3 6 5 3 1 1 5 4 4; 2 4 5 2 5 3 4 4 6; 2 3 3 3 6 4 5 6 5; 3 3 6 2 1 6 3 4 4], [2 6 4 5 6 3 3 4 3; 1 1 4 1 5 2 3 3 2; 3 5 2 2 6 3 1 4 3; 2 6 4 5 1 3 6 5 5; 2 3 6 2 5 6 2 6 3; 4 4 1 6 2 2 6 4 1; 4 3 5 4 1 2 6 1 5; 2 2 1 6 1 4 3 2 1; 3 6 4 2 2 1 4 6 4]]


# starting_configuration_L_11 = [[1 4 3 2 2 4 3 5 2 1 5; 1 3 1 3 1 2 2 1 4 3 1; 1 3 4 5 5 1 3 5 5 2 5; 1 5 6 6 4 2 1 2 6 4 3; 2 2 1 5 2 5 6 6 2 6 5; 2 4 6 5 2 5 4 5 2 3 6; 4 1 6 3 6 6 6 4 5 5 5; 5 4 5 2 3 3 2 2 4 1 6; 4 5 1 6 6 1 6 5 4 2 2; 5 2 2 1 4 1 6 5 2 5 2; 6 3 6 6 2 1 3 5 6 2 4],[4 2 2 1 3 2 5 3 5 3 1; 4 4 6 4 2 6 5 6 1 2 1; 1 1 4 2 3 1 6 2 6 2 6; 3 2 3 6 5 5 4 5 5 4 1; 3 1 2 1 1 3 3 6 3 4 5; 6 2 5 1 6 6 6 3 1 6 3; 5 3 6 5 6 2 1 3 5 3 6; 5 1 4 6 2 5 1 1 1 5 3; 5 2 4 6 4 2 4 2 5 2 5; 4 1 5 4 6 1 4 6 2 2 1; 3 5 2 6 6 1 1 4 5 2 3],[5 2 1 1 1 3 4 5 2 4 3; 6 2 3 4 5 5 2 3 6 6 6; 3 6 5 5 2 1 5 2 4 1 2; 1 2 1 2 3 1 6 4 5 2 5; 6 1 3 4 4 1 2 3 3 4 5; 1 3 5 3 6 3 1 5 1 3 1; 3 6 3 2 1 3 6 5 2 6 1; 4 3 2 3 5 4 2 1 2 5 2; 5 6 1 6 2 6 6 4 5 2 1; 3 1 6 5 5 4 3 5 1 5 4; 6 3 3 1 5 1 4 3 6 5 4],[6 4 1 5 6 3 6 3 2 3 3; 3 5 3 3 6 1 1 6 4 1 1; 1 2 2 5 4 5 3 4 4 1 3; 6 1 6 5 4 5 1 6 6 6 6; 4 3 1 2 4 1 1 5 5 3 6; 4 5 3 3 6 4 6 1 4 4 4; 2 4 6 2 2 1 4 6 3 6 6; 5 4 4 6 1 1 3 4 1 5 4; 4 5 4 4 4 4 3 5 1 6 5; 5 2 4 2 2 5 3 3 3 6 3; 3 2 5 5 5 5 6 5 3 2 6],[6 4 2 2 2 1 3 3 1 1 3; 1 6 4 5 3 3 2 3 6 1 2; 6 4 2 4 3 5 4 2 2 6 2; 3 1 6 5 6 4 6 5 4 2 5; 6 2 4 2 2 3 5 4 1 2 2; 4 6 1 3 3 5 3 3 4 1 2; 3 2 4 1 6 3 4 6 4 6 3; 4 3 5 5 1 1 3 6 1 5 3; 3 3 5 5 6 4 2 1 6 4 1; 1 4 1 6 5 6 3 3 2 1 5; 5 3 2 4 5 2 1 5 2 5 1],[2 3 2 2 1 4 1 2 3 4 6; 6 2 4 4 3 4 1 4 6 4 5; 1 2 3 1 1 1 6 3 1 4 3; 3 1 6 4 6 1 4 3 6 4 3; 3 1 3 4 2 2 6 2 2 1 2; 2 4 6 2 5 3 5 6 3 4 3; 6 2 4 5 4 2 6 1 4 6 4; 4 3 6 5 5 1 3 4 5 4 4; 4 3 6 4 5 2 5 4 5 3 2; 3 1 4 4 5 5 1 1 5 2 6; 2 4 1 4 6 2 6 2 6 4 3]]
# starting_configuration_L_13 = [[6 4 5 1 5 4 2 3 6 6 4 5 5; 3 3 5 3 4 4 4 5 4 1 2 1 3; 6 3 6 4 6 3 4 1 1 5 2 6 3; 3 1 2 3 2 4 5 3 5 6 1 6 5; 5 1 6 3 3 4 6 4 2 1 2 2 2; 4 4 6 3 5 5 1 3 2 4 6 5 1; 4 5 3 4 6 3 1 6 4 5 2 6 3; 5 2 1 3 1 1 1 1 5 2 4 4 6; 6 3 5 3 2 1 5 3 1 3 4 5 6; 3 2 5 4 6 3 2 1 5 2 6 6 1; 6 4 2 5 1 1 1 3 5 3 3 2 2; 5 5 4 6 6 1 1 4 5 5 4 2 1; 2 1 1 1 4 1 2 1 2 4 5 6 4],[5 3 5 4 4 6 5 3 2 6 4 2 5; 5 6 3 1 1 3 2 6 3 6 1 1 5; 1 6 3 4 2 4 5 3 3 6 3 6 6; 4 5 4 6 5 4 6 4 2 2 5 1 3; 3 3 6 1 6 1 1 2 1 4 6 4 5; 1 6 1 5 3 5 3 5 4 4 3 5 6; 3 2 2 4 2 6 4 1 5 6 2 3 3; 4 6 1 4 1 1 4 3 2 1 4 3 1; 2 2 4 2 2 2 2 6 3 5 5 6 6; 2 5 1 5 6 2 2 2 3 6 5 4 3; 3 6 6 4 5 1 2 2 1 2 2 6 4; 1 1 4 4 6 2 2 5 6 2 2 2 2; 3 6 3 3 1 3 4 3 1 4 3 2 1],[3 6 2 5 6 4 5 6 4 4 6 2 3; 1 1 2 4 1 5 3 5 6 1 1 1 5; 5 3 4 3 4 5 2 4 4 3 5 5 4; 6 6 6 4 2 6 4 6 5 6 2 6 2; 2 6 3 5 5 1 6 6 1 6 3 3 6; 1 5 4 6 4 2 2 1 4 1 2 4 6; 2 1 4 4 5 5 1 6 1 5 1 4 2; 4 1 2 3 5 1 6 4 6 3 4 4 3; 1 5 5 2 4 1 1 1 3 3 5 6 6; 6 4 1 1 6 2 1 1 2 5 6 4 1; 5 3 3 3 5 3 1 1 3 5 2 3 3; 6 2 2 1 2 4 2 1 2 1 3 1 3; 4 4 3 2 6 5 3 2 6 2 3 3 5],[6 1 3 2 1 1 4 3 4 2 5 2 3; 2 5 4 2 4 3 4 2 6 2 3 5 2; 3 1 3 1 1 5 1 2 1 3 2 4 4; 4 1 4 6 6 2 4 1 1 5 3 4 3; 3 1 5 4 5 1 3 5 4 5 1 5 2; 1 6 6 2 6 3 5 5 3 3 4 5 5; 4 1 1 5 6 3 4 3 6 2 1 2 1; 1 5 3 6 3 4 1 6 1 5 1 6 4; 6 5 2 6 5 5 2 3 6 4 5 6 5; 5 5 3 5 3 6 4 6 4 5 3 5 6; 5 2 3 4 6 3 1 6 6 6 5 1 4; 6 3 4 5 4 6 3 4 5 2 2 6 6; 4 5 1 2 1 2 6 3 4 1 3 2 2],[5 4 2 3 1 5 2 6 3 3 2 4 5; 1 1 2 6 5 6 2 1 3 4 1 3 5; 2 2 6 4 6 2 5 5 4 4 2 6 1; 6 2 2 4 5 3 3 2 6 4 6 6 1; 1 1 4 4 6 5 2 5 4 1 3 4 2; 4 5 1 6 1 6 2 4 5 3 2 6 2; 2 5 4 4 4 1 6 1 3 5 4 4 2; 4 5 3 1 2 3 2 4 1 1 4 3 3; 5 4 3 1 4 3 3 5 1 3 4 3 4; 1 2 5 6 3 6 6 4 6 6 6 2 2; 1 5 3 1 2 1 5 2 4 2 6 2 1; 6 2 1 5 4 6 1 3 6 4 6 5 2; 5 2 5 3 4 4 2 2 1 5 1 2 4],[3 5 5 1 6 6 3 4 2 6 3 1 3; 4 5 4 3 5 6 5 1 4 6 6 4 2; 4 6 5 3 6 5 2 3 6 2 2 4 3; 1 4 1 4 3 5 6 2 3 5 6 5 2; 5 2 3 2 3 2 6 3 3 5 1 5 2; 6 2 1 3 5 4 5 3 2 6 4 6 5; 4 2 5 4 3 2 2 5 3 6 3 6 5; 2 3 3 5 2 6 6 2 5 3 5 6 3; 2 1 6 2 3 5 4 2 3 1 1 5 3; 4 3 4 5 5 4 1 5 6 6 2 5 6; 3 3 1 1 4 1 5 4 1 6 2 4 3; 2 6 1 2 1 5 5 3 4 4 4 1 2; 5 6 1 4 6 1 1 3 5 6 2 1 1]]
# starting_configuration_L_15 = [[3 5 1 1 3 3 5 3 6 4 4 6 6 5 3; 2 1 4 3 6 1 4 5 2 4 1 5 6 3 3; 5 2 1 6 2 4 1 6 5 3 3 6 4 2 2; 2 5 5 1 5 3 2 5 5 1 5 1 6 5 2; 2 6 6 6 5 2 2 6 4 5 6 3 1 4 4; 3 2 5 3 3 1 5 5 1 4 2 2 1 6 3; 6 6 3 6 4 6 2 3 3 5 5 2 3 2 6; 2 6 5 5 5 2 5 6 2 2 2 5 1 6 6; 4 6 2 4 5 1 5 6 6 6 4 1 5 1 5; 5 6 1 6 1 2 5 6 3 6 6 1 1 2 3; 1 5 5 1 2 3 5 4 3 6 4 1 4 2 4; 3 2 1 2 6 2 2 5 4 2 3 3 1 6 6; 4 5 1 3 2 1 2 2 4 1 6 4 4 2 1; 1 6 3 1 1 4 3 3 3 6 4 2 3 4 3; 2 1 2 1 3 5 1 5 3 6 4 3 2 1 1],[1 3 6 3 3 3 4 5 3 2 5 1 6 6 3; 5 1 4 6 1 5 5 6 6 1 6 4 4 1 3; 6 1 1 2 5 6 4 1 3 2 4 6 5 2 5; 1 3 4 1 4 3 1 3 5 6 3 2 4 1 1; 2 6 6 5 3 5 2 1 4 2 6 4 6 3 6; 4 6 1 6 2 1 3 5 5 4 4 2 6 6 5; 6 3 1 6 1 3 4 4 3 3 5 1 6 1 6; 6 5 5 3 5 6 2 1 3 2 2 3 1 5 4; 1 1 2 4 3 3 1 5 1 1 5 6 4 5 2; 2 4 6 1 5 4 3 2 4 5 2 4 6 5 5; 1 6 5 4 1 3 6 1 3 5 3 3 1 1 5; 2 3 4 3 2 1 6 6 5 3 3 5 4 5 1; 1 3 2 2 4 2 2 6 2 6 5 3 5 5 3; 6 3 4 3 3 3 5 5 4 5 3 2 6 1 5; 1 3 4 3 4 4 1 6 4 1 6 3 2 6 1],[4 2 5 1 6 5 1 1 2 1 2 3 3 1 5; 1 2 1 5 1 2 3 4 2 3 3 1 5 2 3; 3 1 6 6 1 1 5 3 2 4 6 3 1 4 4; 5 5 5 2 4 1 3 6 5 5 5 4 4 3 1; 1 4 5 1 2 4 6 4 6 4 3 2 5 2 1; 2 3 5 1 3 1 3 1 4 6 4 6 3 3 5; 3 2 4 6 5 6 2 2 5 4 4 5 1 3 3; 6 6 6 4 1 2 6 2 5 1 1 2 2 2 1; 3 2 6 6 2 1 6 4 6 2 3 5 1 2 5; 5 5 2 5 4 3 2 3 5 5 2 1 6 4 5; 5 5 5 3 5 5 5 1 3 6 4 6 6 6 2; 2 3 6 4 5 4 3 5 1 4 5 2 3 6 2; 3 3 2 5 1 5 3 4 5 6 1 1 2 6 1; 6 5 5 2 1 6 4 1 5 4 4 1 4 5 3; 2 6 2 4 4 5 1 3 6 4 5 1 1 3 2],[6 1 6 6 4 2 3 2 6 4 6 4 6 1 1; 3 1 2 2 4 1 6 1 2 5 5 6 5 1 3; 6 2 4 1 2 3 2 2 4 3 6 5 2 6 6; 2 5 5 1 6 5 5 4 3 2 4 5 3 4 2; 5 1 4 1 3 2 3 3 6 2 6 4 6 3 3; 2 2 6 4 3 6 1 6 1 2 2 4 4 3 1; 5 4 4 6 1 5 1 1 4 4 2 5 3 6 6; 1 2 3 5 2 3 6 6 5 4 2 1 5 3 3; 6 4 2 3 3 6 2 2 6 6 3 3 6 3 4; 6 4 3 5 1 5 6 4 2 3 6 4 2 6 3; 5 5 6 1 1 4 1 4 4 4 2 6 4 2 3; 3 2 4 2 1 1 3 2 2 1 2 5 2 4 4; 3 4 2 4 6 6 2 1 4 6 5 3 2 6 5; 4 6 2 3 6 6 3 1 3 6 1 6 4 1 6; 6 5 4 4 2 6 4 2 1 6 6 2 1 3 4],[3 1 2 1 6 4 3 2 6 4 3 3 4 5 4; 4 3 2 3 2 3 6 3 3 2 3 2 1 3 2; 4 4 3 6 6 1 1 6 3 3 3 2 4 5 5; 1 5 5 2 4 3 6 1 6 4 1 6 6 4 6; 2 3 3 4 1 3 4 5 3 4 3 6 2 2 4; 6 4 5 5 1 1 5 6 2 1 1 4 4 2 3; 5 1 1 5 1 4 3 5 5 4 2 6 5 4 6; 4 6 2 5 4 6 5 2 5 4 4 2 4 6 3; 5 4 2 4 2 1 5 4 4 4 4 5 1 6 5; 4 4 6 1 5 2 1 3 3 3 2 5 5 5 4; 6 3 2 1 5 1 4 5 5 2 1 6 3 1 5; 2 4 4 3 3 3 4 5 1 6 1 6 2 6 1; 3 2 3 6 3 1 1 3 4 3 5 2 3 5 2; 3 1 5 6 1 6 1 3 5 3 1 2 2 3 6; 2 6 2 5 5 5 5 6 2 1 1 1 6 6 5],[4 2 2 1 5 1 2 2 4 4 4 3 5 6 3; 1 2 4 1 3 4 6 6 4 6 2 4 6 6 1; 4 4 1 3 2 4 5 4 2 4 5 1 2 2 1; 5 4 5 6 5 2 3 2 3 3 6 5 2 1 4; 4 2 4 1 3 5 4 3 5 3 4 2 6 4 2; 5 2 1 3 5 2 4 5 3 3 2 1 4 4 5; 5 6 1 2 4 6 2 5 1 5 3 3 6 1 1; 3 4 1 2 5 3 2 6 2 3 4 5 6 1 6; 3 4 1 5 5 4 4 1 2 6 4 3 4 4 3; 2 2 4 5 4 4 2 1 1 1 6 4 6 2 6; 4 4 2 4 3 4 5 5 6 3 3 3 2 1 5; 1 1 5 5 2 2 4 6 3 5 6 5 4 4 4; 6 2 6 5 6 5 1 5 1 2 6 1 4 5 2; 6 2 4 3 1 6 3 2 2 5 2 5 1 1 4; 3 4 2 4 5 1 3 1 6 5 3 4 2 1 5]]

# L_values = [11,13,15]
# starting_configurations = [starting_configuration_L_11, starting_configuration_L_13, starting_configuration_L_15]

L_values = [5,7]
glassy_configurations = [glassy_configuration_L_5, glassy_configuration_L_7]
energy_caps = [-100.0,-230.0]


@sync @distributed for L_value in L_values
    capped_anneal_experiment("capped_anneal_L=$L_value", glassy_configurations[Int((L_value-3)/2)], energy_caps[Int((L_value-3)/2)], 0.0, 1.0, 0.3, Int(2e8); infinite_temperature_steps=50, verbose_annealing=true)
end