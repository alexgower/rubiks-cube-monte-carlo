using Distributed
addprocs(3)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/relaxed_anneal_experiment.jl")


starting_configuration_L_5 = [[2 3 3 2 6; 5 3 6 4 4; 1 2 5 1 5; 2 2 6 6 2; 2 1 2 3 2],[1 2 5 5 3; 6 6 1 2 3; 2 5 5 4 3; 4 1 6 5 4; 4 6 6 2 4], [3 4 5 4 6; 2 5 4 4 1; 5 3 2 5 4; 4 6 1 2 5; 2 5 5 4 6], [3 3 1 6 3; 2 3 1 6 3; 6 6 4 3 2; 1 1 2 3 2; 6 3 1 5 3], [2 3 4 4 1; 3 3 1 3 6; 3 4 4 5 1; 4 5 5 4 1; 6 1 1 6 6], [4 1 2 3 5; 2 1 6 6 5; 4 1 5 5 1; 2 1 5 1 6; 5 4 4 3 6]]
starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 4 1 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 6 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 2; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 6 5 2 2; 3 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 4 3 1 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 2 1 6 3; 6 3 6 6 5 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
starting_configuration_L_9 = [[2 4 2 6 4 5 2 6 2; 6 6 3 5 1 3 4 3 4; 5 6 6 3 1 5 1 4 4; 4 1 5 6 6 3 1 4 1; 2 6 4 4 6 2 4 5 5; 5 5 1 2 2 2 1 4 3; 3 4 1 1 5 3 3 3 6; 3 4 2 3 2 4 3 4 2; 6 6 5 1 1 3 5 4 4], [6 6 6 3 2 6 1 2 4; 5 1 5 1 2 6 6 2 4; 4 5 2 3 3 1 5 1 2; 4 4 3 4 6 1 6 5 1; 5 5 4 6 6 3 4 3 6; 5 5 4 6 2 3 4 1 4; 6 6 6 2 4 6 1 6 6; 5 5 1 5 5 1 5 4 1; 4 3 5 5 4 1 4 2 5], [5 2 5 2 6 3 1 2 2; 5 3 3 4 3 5 5 4 1; 3 5 4 2 3 1 6 5 2; 4 5 3 2 5 4 2 2 3; 5 4 5 5 1 3 2 6 1; 6 4 3 6 5 1 3 1 6; 5 4 6 3 1 5 3 2 1; 6 5 1 4 3 1 3 6 3; 2 4 5 6 2 1 1 3 5], [3 2 3 5 4 6 2 6 3; 4 6 3 4 3 2 3 1 5; 4 1 2 1 2 4 2 5 3; 6 6 6 1 1 2 6 6 5; 2 2 6 6 6 2 4 1 3; 1 1 2 1 4 3 2 1 6; 3 2 1 1 4 1 4 3 2; 5 3 2 3 1 3 2 1 3; 5 5 1 6 1 5 5 1 2], [5 3 4 6 4 2 5 5 6; 4 2 5 5 1 1 5 2 1; 6 1 2 4 6 2 1 1 4; 6 3 1 5 1 6 2 5 6; 1 2 2 3 3 2 5 5 4; 3 6 5 3 1 1 5 4 4; 2 4 5 2 5 3 4 4 6; 2 3 3 3 6 4 5 6 5; 3 3 6 2 1 6 3 4 4], [2 6 4 5 6 3 3 4 3; 1 1 4 1 5 2 3 3 2; 3 5 2 2 6 3 1 4 3; 2 6 4 5 1 3 6 5 5; 2 3 6 2 5 6 2 6 3; 4 4 1 6 2 2 6 4 1; 4 3 5 4 1 2 6 1 5; 2 2 1 6 1 4 3 2 1; 3 6 4 2 2 1 4 6 4]]

new_starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 1 4 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 2 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 6; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 3 5 2 2; 6 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 1 3 1 4 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 5 1 6 3; 6 3 6 6 2 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 1 2 4 6 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
new_new_starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 1 4 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5], [1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 2 1 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 6; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 3 1 5 2 2; 6 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 1 3 4 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 5 1 6 3; 6 3 6 6 2 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]


# starting_configuration_L_11 = [[1 4 3 2 2 4 3 5 2 1 5; 1 3 1 3 1 2 2 1 4 3 1; 1 3 4 5 5 1 3 5 5 2 5; 1 5 6 6 4 2 1 2 6 4 3; 2 2 1 5 2 5 6 6 2 6 5; 2 4 6 5 2 5 4 5 2 3 6; 4 1 6 3 6 6 6 4 5 5 5; 5 4 5 2 3 3 2 2 4 1 6; 4 5 1 6 6 1 6 5 4 2 2; 5 2 2 1 4 1 6 5 2 5 2; 6 3 6 6 2 1 3 5 6 2 4],[4 2 2 1 3 2 5 3 5 3 1; 4 4 6 4 2 6 5 6 1 2 1; 1 1 4 2 3 1 6 2 6 2 6; 3 2 3 6 5 5 4 5 5 4 1; 3 1 2 1 1 3 3 6 3 4 5; 6 2 5 1 6 6 6 3 1 6 3; 5 3 6 5 6 2 1 3 5 3 6; 5 1 4 6 2 5 1 1 1 5 3; 5 2 4 6 4 2 4 2 5 2 5; 4 1 5 4 6 1 4 6 2 2 1; 3 5 2 6 6 1 1 4 5 2 3],[5 2 1 1 1 3 4 5 2 4 3; 6 2 3 4 5 5 2 3 6 6 6; 3 6 5 5 2 1 5 2 4 1 2; 1 2 1 2 3 1 6 4 5 2 5; 6 1 3 4 4 1 2 3 3 4 5; 1 3 5 3 6 3 1 5 1 3 1; 3 6 3 2 1 3 6 5 2 6 1; 4 3 2 3 5 4 2 1 2 5 2; 5 6 1 6 2 6 6 4 5 2 1; 3 1 6 5 5 4 3 5 1 5 4; 6 3 3 1 5 1 4 3 6 5 4],[6 4 1 5 6 3 6 3 2 3 3; 3 5 3 3 6 1 1 6 4 1 1; 1 2 2 5 4 5 3 4 4 1 3; 6 1 6 5 4 5 1 6 6 6 6; 4 3 1 2 4 1 1 5 5 3 6; 4 5 3 3 6 4 6 1 4 4 4; 2 4 6 2 2 1 4 6 3 6 6; 5 4 4 6 1 1 3 4 1 5 4; 4 5 4 4 4 4 3 5 1 6 5; 5 2 4 2 2 5 3 3 3 6 3; 3 2 5 5 5 5 6 5 3 2 6],[6 4 2 2 2 1 3 3 1 1 3; 1 6 4 5 3 3 2 3 6 1 2; 6 4 2 4 3 5 4 2 2 6 2; 3 1 6 5 6 4 6 5 4 2 5; 6 2 4 2 2 3 5 4 1 2 2; 4 6 1 3 3 5 3 3 4 1 2; 3 2 4 1 6 3 4 6 4 6 3; 4 3 5 5 1 1 3 6 1 5 3; 3 3 5 5 6 4 2 1 6 4 1; 1 4 1 6 5 6 3 3 2 1 5; 5 3 2 4 5 2 1 5 2 5 1],[2 3 2 2 1 4 1 2 3 4 6; 6 2 4 4 3 4 1 4 6 4 5; 1 2 3 1 1 1 6 3 1 4 3; 3 1 6 4 6 1 4 3 6 4 3; 3 1 3 4 2 2 6 2 2 1 2; 2 4 6 2 5 3 5 6 3 4 3; 6 2 4 5 4 2 6 1 4 6 4; 4 3 6 5 5 1 3 4 5 4 4; 4 3 6 4 5 2 5 4 5 3 2; 3 1 4 4 5 5 1 1 5 2 6; 2 4 1 4 6 2 6 2 6 4 3]]
# starting_configuration_L_13 = [[6 4 5 1 5 4 2 3 6 6 4 5 5; 3 3 5 3 4 4 4 5 4 1 2 1 3; 6 3 6 4 6 3 4 1 1 5 2 6 3; 3 1 2 3 2 4 5 3 5 6 1 6 5; 5 1 6 3 3 4 6 4 2 1 2 2 2; 4 4 6 3 5 5 1 3 2 4 6 5 1; 4 5 3 4 6 3 1 6 4 5 2 6 3; 5 2 1 3 1 1 1 1 5 2 4 4 6; 6 3 5 3 2 1 5 3 1 3 4 5 6; 3 2 5 4 6 3 2 1 5 2 6 6 1; 6 4 2 5 1 1 1 3 5 3 3 2 2; 5 5 4 6 6 1 1 4 5 5 4 2 1; 2 1 1 1 4 1 2 1 2 4 5 6 4],[5 3 5 4 4 6 5 3 2 6 4 2 5; 5 6 3 1 1 3 2 6 3 6 1 1 5; 1 6 3 4 2 4 5 3 3 6 3 6 6; 4 5 4 6 5 4 6 4 2 2 5 1 3; 3 3 6 1 6 1 1 2 1 4 6 4 5; 1 6 1 5 3 5 3 5 4 4 3 5 6; 3 2 2 4 2 6 4 1 5 6 2 3 3; 4 6 1 4 1 1 4 3 2 1 4 3 1; 2 2 4 2 2 2 2 6 3 5 5 6 6; 2 5 1 5 6 2 2 2 3 6 5 4 3; 3 6 6 4 5 1 2 2 1 2 2 6 4; 1 1 4 4 6 2 2 5 6 2 2 2 2; 3 6 3 3 1 3 4 3 1 4 3 2 1],[3 6 2 5 6 4 5 6 4 4 6 2 3; 1 1 2 4 1 5 3 5 6 1 1 1 5; 5 3 4 3 4 5 2 4 4 3 5 5 4; 6 6 6 4 2 6 4 6 5 6 2 6 2; 2 6 3 5 5 1 6 6 1 6 3 3 6; 1 5 4 6 4 2 2 1 4 1 2 4 6; 2 1 4 4 5 5 1 6 1 5 1 4 2; 4 1 2 3 5 1 6 4 6 3 4 4 3; 1 5 5 2 4 1 1 1 3 3 5 6 6; 6 4 1 1 6 2 1 1 2 5 6 4 1; 5 3 3 3 5 3 1 1 3 5 2 3 3; 6 2 2 1 2 4 2 1 2 1 3 1 3; 4 4 3 2 6 5 3 2 6 2 3 3 5],[6 1 3 2 1 1 4 3 4 2 5 2 3; 2 5 4 2 4 3 4 2 6 2 3 5 2; 3 1 3 1 1 5 1 2 1 3 2 4 4; 4 1 4 6 6 2 4 1 1 5 3 4 3; 3 1 5 4 5 1 3 5 4 5 1 5 2; 1 6 6 2 6 3 5 5 3 3 4 5 5; 4 1 1 5 6 3 4 3 6 2 1 2 1; 1 5 3 6 3 4 1 6 1 5 1 6 4; 6 5 2 6 5 5 2 3 6 4 5 6 5; 5 5 3 5 3 6 4 6 4 5 3 5 6; 5 2 3 4 6 3 1 6 6 6 5 1 4; 6 3 4 5 4 6 3 4 5 2 2 6 6; 4 5 1 2 1 2 6 3 4 1 3 2 2],[5 4 2 3 1 5 2 6 3 3 2 4 5; 1 1 2 6 5 6 2 1 3 4 1 3 5; 2 2 6 4 6 2 5 5 4 4 2 6 1; 6 2 2 4 5 3 3 2 6 4 6 6 1; 1 1 4 4 6 5 2 5 4 1 3 4 2; 4 5 1 6 1 6 2 4 5 3 2 6 2; 2 5 4 4 4 1 6 1 3 5 4 4 2; 4 5 3 1 2 3 2 4 1 1 4 3 3; 5 4 3 1 4 3 3 5 1 3 4 3 4; 1 2 5 6 3 6 6 4 6 6 6 2 2; 1 5 3 1 2 1 5 2 4 2 6 2 1; 6 2 1 5 4 6 1 3 6 4 6 5 2; 5 2 5 3 4 4 2 2 1 5 1 2 4],[3 5 5 1 6 6 3 4 2 6 3 1 3; 4 5 4 3 5 6 5 1 4 6 6 4 2; 4 6 5 3 6 5 2 3 6 2 2 4 3; 1 4 1 4 3 5 6 2 3 5 6 5 2; 5 2 3 2 3 2 6 3 3 5 1 5 2; 6 2 1 3 5 4 5 3 2 6 4 6 5; 4 2 5 4 3 2 2 5 3 6 3 6 5; 2 3 3 5 2 6 6 2 5 3 5 6 3; 2 1 6 2 3 5 4 2 3 1 1 5 3; 4 3 4 5 5 4 1 5 6 6 2 5 6; 3 3 1 1 4 1 5 4 1 6 2 4 3; 2 6 1 2 1 5 5 3 4 4 4 1 2; 5 6 1 4 6 1 1 3 5 6 2 1 1]]
# starting_configuration_L_15 = [[3 5 1 1 3 3 5 3 6 4 4 6 6 5 3; 2 1 4 3 6 1 4 5 2 4 1 5 6 3 3; 5 2 1 6 2 4 1 6 5 3 3 6 4 2 2; 2 5 5 1 5 3 2 5 5 1 5 1 6 5 2; 2 6 6 6 5 2 2 6 4 5 6 3 1 4 4; 3 2 5 3 3 1 5 5 1 4 2 2 1 6 3; 6 6 3 6 4 6 2 3 3 5 5 2 3 2 6; 2 6 5 5 5 2 5 6 2 2 2 5 1 6 6; 4 6 2 4 5 1 5 6 6 6 4 1 5 1 5; 5 6 1 6 1 2 5 6 3 6 6 1 1 2 3; 1 5 5 1 2 3 5 4 3 6 4 1 4 2 4; 3 2 1 2 6 2 2 5 4 2 3 3 1 6 6; 4 5 1 3 2 1 2 2 4 1 6 4 4 2 1; 1 6 3 1 1 4 3 3 3 6 4 2 3 4 3; 2 1 2 1 3 5 1 5 3 6 4 3 2 1 1],[1 3 6 3 3 3 4 5 3 2 5 1 6 6 3; 5 1 4 6 1 5 5 6 6 1 6 4 4 1 3; 6 1 1 2 5 6 4 1 3 2 4 6 5 2 5; 1 3 4 1 4 3 1 3 5 6 3 2 4 1 1; 2 6 6 5 3 5 2 1 4 2 6 4 6 3 6; 4 6 1 6 2 1 3 5 5 4 4 2 6 6 5; 6 3 1 6 1 3 4 4 3 3 5 1 6 1 6; 6 5 5 3 5 6 2 1 3 2 2 3 1 5 4; 1 1 2 4 3 3 1 5 1 1 5 6 4 5 2; 2 4 6 1 5 4 3 2 4 5 2 4 6 5 5; 1 6 5 4 1 3 6 1 3 5 3 3 1 1 5; 2 3 4 3 2 1 6 6 5 3 3 5 4 5 1; 1 3 2 2 4 2 2 6 2 6 5 3 5 5 3; 6 3 4 3 3 3 5 5 4 5 3 2 6 1 5; 1 3 4 3 4 4 1 6 4 1 6 3 2 6 1],[4 2 5 1 6 5 1 1 2 1 2 3 3 1 5; 1 2 1 5 1 2 3 4 2 3 3 1 5 2 3; 3 1 6 6 1 1 5 3 2 4 6 3 1 4 4; 5 5 5 2 4 1 3 6 5 5 5 4 4 3 1; 1 4 5 1 2 4 6 4 6 4 3 2 5 2 1; 2 3 5 1 3 1 3 1 4 6 4 6 3 3 5; 3 2 4 6 5 6 2 2 5 4 4 5 1 3 3; 6 6 6 4 1 2 6 2 5 1 1 2 2 2 1; 3 2 6 6 2 1 6 4 6 2 3 5 1 2 5; 5 5 2 5 4 3 2 3 5 5 2 1 6 4 5; 5 5 5 3 5 5 5 1 3 6 4 6 6 6 2; 2 3 6 4 5 4 3 5 1 4 5 2 3 6 2; 3 3 2 5 1 5 3 4 5 6 1 1 2 6 1; 6 5 5 2 1 6 4 1 5 4 4 1 4 5 3; 2 6 2 4 4 5 1 3 6 4 5 1 1 3 2],[6 1 6 6 4 2 3 2 6 4 6 4 6 1 1; 3 1 2 2 4 1 6 1 2 5 5 6 5 1 3; 6 2 4 1 2 3 2 2 4 3 6 5 2 6 6; 2 5 5 1 6 5 5 4 3 2 4 5 3 4 2; 5 1 4 1 3 2 3 3 6 2 6 4 6 3 3; 2 2 6 4 3 6 1 6 1 2 2 4 4 3 1; 5 4 4 6 1 5 1 1 4 4 2 5 3 6 6; 1 2 3 5 2 3 6 6 5 4 2 1 5 3 3; 6 4 2 3 3 6 2 2 6 6 3 3 6 3 4; 6 4 3 5 1 5 6 4 2 3 6 4 2 6 3; 5 5 6 1 1 4 1 4 4 4 2 6 4 2 3; 3 2 4 2 1 1 3 2 2 1 2 5 2 4 4; 3 4 2 4 6 6 2 1 4 6 5 3 2 6 5; 4 6 2 3 6 6 3 1 3 6 1 6 4 1 6; 6 5 4 4 2 6 4 2 1 6 6 2 1 3 4],[3 1 2 1 6 4 3 2 6 4 3 3 4 5 4; 4 3 2 3 2 3 6 3 3 2 3 2 1 3 2; 4 4 3 6 6 1 1 6 3 3 3 2 4 5 5; 1 5 5 2 4 3 6 1 6 4 1 6 6 4 6; 2 3 3 4 1 3 4 5 3 4 3 6 2 2 4; 6 4 5 5 1 1 5 6 2 1 1 4 4 2 3; 5 1 1 5 1 4 3 5 5 4 2 6 5 4 6; 4 6 2 5 4 6 5 2 5 4 4 2 4 6 3; 5 4 2 4 2 1 5 4 4 4 4 5 1 6 5; 4 4 6 1 5 2 1 3 3 3 2 5 5 5 4; 6 3 2 1 5 1 4 5 5 2 1 6 3 1 5; 2 4 4 3 3 3 4 5 1 6 1 6 2 6 1; 3 2 3 6 3 1 1 3 4 3 5 2 3 5 2; 3 1 5 6 1 6 1 3 5 3 1 2 2 3 6; 2 6 2 5 5 5 5 6 2 1 1 1 6 6 5],[4 2 2 1 5 1 2 2 4 4 4 3 5 6 3; 1 2 4 1 3 4 6 6 4 6 2 4 6 6 1; 4 4 1 3 2 4 5 4 2 4 5 1 2 2 1; 5 4 5 6 5 2 3 2 3 3 6 5 2 1 4; 4 2 4 1 3 5 4 3 5 3 4 2 6 4 2; 5 2 1 3 5 2 4 5 3 3 2 1 4 4 5; 5 6 1 2 4 6 2 5 1 5 3 3 6 1 1; 3 4 1 2 5 3 2 6 2 3 4 5 6 1 6; 3 4 1 5 5 4 4 1 2 6 4 3 4 4 3; 2 2 4 5 4 4 2 1 1 1 6 4 6 2 6; 4 4 2 4 3 4 5 5 6 3 3 3 2 1 5; 1 1 5 5 2 2 4 6 3 5 6 5 4 4 4; 6 2 6 5 6 5 1 5 1 2 6 1 4 5 2; 6 2 4 3 1 6 3 2 2 5 2 5 1 1 4; 3 4 2 4 5 1 3 1 6 5 3 4 2 1 5]]

# L_values = [11,13,15]
# starting_configurations = [starting_configuration_L_11, starting_configuration_L_13, starting_configuration_L_15]

# L_values = [5,7,9]
# starting_configurations = [starting_configuration_L_5, starting_configuration_L_7, starting_configuration_L_9]

# TODO CHECK IF WANT NEW OR OLD STARTING CONFIGURATION FOR L=7
L_values = [7]
starting_configurations = [starting_configuration_L_5, new_new_starting_configuration_L_7, starting_configuration_L_9]


@sync @distributed for L_value in L_values
    # TODO ADD ps 1.0 back
    relaxed_anneal_experiment("relaxed_anneal_new_new_L=$L_value", L_value, [0.0], 10.0,10.0,0.1,120; inherent_disorder=true, relaxation_iterations=100000, initial_cube_configuration=starting_configurations[Int((L_value-3)/2)])
end