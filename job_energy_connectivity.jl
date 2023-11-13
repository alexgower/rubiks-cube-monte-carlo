number_of_processors = 50
L_values = [9,10]

number_of_L_values = length(L_values)
number_of_processors_per_L_value = Int(floor(number_of_processors/number_of_L_values))


using Distributed
addprocs(50)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/neighbour_initial_and_final_energies_distribution_experiment.jl")

starting_configuration_L_5 = [[2 3 3 2 6; 5 3 6 4 4; 1 2 5 1 5; 2 2 6 6 2; 2 1 2 3 2],[1 2 5 5 3; 6 6 1 2 3; 2 5 5 4 3; 4 1 6 5 4; 4 6 6 2 4], [3 4 5 4 6; 2 5 4 4 1; 5 3 2 5 4; 4 6 1 2 5; 2 5 5 4 6], [3 3 1 6 3; 2 3 1 6 3; 6 6 4 3 2; 1 1 2 3 2; 6 3 1 5 3], [2 3 4 4 1; 3 3 1 3 6; 3 4 4 5 1; 4 5 5 4 1; 6 1 1 6 6], [4 1 2 3 5; 2 1 6 6 5; 4 1 5 5 1; 2 1 5 1 6; 5 4 4 3 6]]
starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 4 1 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 6 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 2; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 6 5 2 2; 3 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 4 3 1 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 2 1 6 3; 6 3 6 6 5 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
starting_configuration_L_9 = [[2 4 2 6 4 5 2 6 2; 6 6 3 5 1 3 4 3 4; 5 6 6 3 1 5 1 4 4; 4 1 5 6 6 3 1 4 1; 2 6 4 4 6 2 4 5 5; 5 5 1 2 2 2 1 4 3; 3 4 1 1 5 3 3 3 6; 3 4 2 3 2 4 3 4 2; 6 6 5 1 1 3 5 4 4], [6 6 6 3 2 6 1 2 4; 5 1 5 1 2 6 6 2 4; 4 5 2 3 3 1 5 1 2; 4 4 3 4 6 1 6 5 1; 5 5 4 6 6 3 4 3 6; 5 5 4 6 2 3 4 1 4; 6 6 6 2 4 6 1 6 6; 5 5 1 5 5 1 5 4 1; 4 3 5 5 4 1 4 2 5], [5 2 5 2 6 3 1 2 2; 5 3 3 4 3 5 5 4 1; 3 5 4 2 3 1 6 5 2; 4 5 3 2 5 4 2 2 3; 5 4 5 5 1 3 2 6 1; 6 4 3 6 5 1 3 1 6; 5 4 6 3 1 5 3 2 1; 6 5 1 4 3 1 3 6 3; 2 4 5 6 2 1 1 3 5], [3 2 3 5 4 6 2 6 3; 4 6 3 4 3 2 3 1 5; 4 1 2 1 2 4 2 5 3; 6 6 6 1 1 2 6 6 5; 2 2 6 6 6 2 4 1 3; 1 1 2 1 4 3 2 1 6; 3 2 1 1 4 1 4 3 2; 5 3 2 3 1 3 2 1 3; 5 5 1 6 1 5 5 1 2], [5 3 4 6 4 2 5 5 6; 4 2 5 5 1 1 5 2 1; 6 1 2 4 6 2 1 1 4; 6 3 1 5 1 6 2 5 6; 1 2 2 3 3 2 5 5 4; 3 6 5 3 1 1 5 4 4; 2 4 5 2 5 3 4 4 6; 2 3 3 3 6 4 5 6 5; 3 3 6 2 1 6 3 4 4], [2 6 4 5 6 3 3 4 3; 1 1 4 1 5 2 3 3 2; 3 5 2 2 6 3 1 4 3; 2 6 4 5 1 3 6 5 5; 2 3 6 2 5 6 2 6 3; 4 4 1 6 2 2 6 4 1; 4 3 5 4 1 2 6 1 5; 2 2 1 6 1 4 3 2 1; 3 6 4 2 2 1 4 6 4]]

L_values = [5,7,9]
starting_configurations = [starting_configuration_L_5, starting_configuration_L_7, starting_configuration_L_9]


@sync @distributed for L_value in L_values
    neighbour_initial_and_final_energies_distribution_experiment("L=$(L_value)_inherent_disorder_E0_E1_slice", L_value, 1.0, 100, [10.0,5.0,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]; relaxation_iterations=10000, collecting_swap_move_neighbours=false, neighbours_per_configuration_sample_size=0, average_sample_size_per_temperature=1000, inherent_disorder=true, neighbour_moves_away=1, parallel_anneals=number_of_processors_per_L_value, initial_cube_configuration=starting_configurations[Int((L_value-3)/2)])
end



