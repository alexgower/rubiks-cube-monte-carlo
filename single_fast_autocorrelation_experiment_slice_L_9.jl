using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("code/experiments/fast_autocorrelation_experiment.jl")



# Command-line arguments
sample_temperature = parse(Float64, ARGS[1])
autocorrelation_window_length = parse(Int, ARGS[2])
trial_number = parse(Int, ARGS[3])

lag_limit = autocorrelation_window_length - 100

# Print out the command-line arguments
println("Sample Temperature: ", sample_temperature)
println("Autocorrelation Window Length: ", autocorrelation_window_length)

println("Running autocorrelation experiment")

# Not using starting configuration now as doing disorder average
# starting_configuration_L_5 = [[2 3 3 2 6; 5 3 6 4 4; 1 2 5 1 5; 2 2 6 6 2; 2 1 2 3 2],[1 2 5 5 3; 6 6 1 2 3; 2 5 5 4 3; 4 1 6 5 4; 4 6 6 2 4], [3 4 5 4 6; 2 5 4 4 1; 5 3 2 5 4; 4 6 1 2 5; 2 5 5 4 6], [3 3 1 6 3; 2 3 1 6 3; 6 6 4 3 2; 1 1 2 3 2; 6 3 1 5 3], [2 3 4 4 1; 3 3 1 3 6; 3 4 4 5 1; 4 5 5 4 1; 6 1 1 6 6], [4 1 2 3 5; 2 1 6 6 5; 4 1 5 5 1; 2 1 5 1 6; 5 4 4 3 6]]
# starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 4 1 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 6 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 2; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 6 5 2 2; 3 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 4 3 1 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 2 1 6 3; 6 3 6 6 5 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
# starting_configuration_L_9 = [[2 4 2 6 4 5 2 6 2; 6 6 3 5 1 3 4 3 4; 5 6 6 3 1 5 1 4 4; 4 1 5 6 6 3 1 4 1; 2 6 4 4 6 2 4 5 5; 5 5 1 2 2 2 1 4 3; 3 4 1 1 5 3 3 3 6; 3 4 2 3 2 4 3 4 2; 6 6 5 1 1 3 5 4 4], [6 6 6 3 2 6 1 2 4; 5 1 5 1 2 6 6 2 4; 4 5 2 3 3 1 5 1 2; 4 4 3 4 6 1 6 5 1; 5 5 4 6 6 3 4 3 6; 5 5 4 6 2 3 4 1 4; 6 6 6 2 4 6 1 6 6; 5 5 1 5 5 1 5 4 1; 4 3 5 5 4 1 4 2 5], [5 2 5 2 6 3 1 2 2; 5 3 3 4 3 5 5 4 1; 3 5 4 2 3 1 6 5 2; 4 5 3 2 5 4 2 2 3; 5 4 5 5 1 3 2 6 1; 6 4 3 6 5 1 3 1 6; 5 4 6 3 1 5 3 2 1; 6 5 1 4 3 1 3 6 3; 2 4 5 6 2 1 1 3 5], [3 2 3 5 4 6 2 6 3; 4 6 3 4 3 2 3 1 5; 4 1 2 1 2 4 2 5 3; 6 6 6 1 1 2 6 6 5; 2 2 6 6 6 2 4 1 3; 1 1 2 1 4 3 2 1 6; 3 2 1 1 4 1 4 3 2; 5 3 2 3 1 3 2 1 3; 5 5 1 6 1 5 5 1 2], [5 3 4 6 4 2 5 5 6; 4 2 5 5 1 1 5 2 1; 6 1 2 4 6 2 1 1 4; 6 3 1 5 1 6 2 5 6; 1 2 2 3 3 2 5 5 4; 3 6 5 3 1 1 5 4 4; 2 4 5 2 5 3 4 4 6; 2 3 3 3 6 4 5 6 5; 3 3 6 2 1 6 3 4 4], [2 6 4 5 6 3 3 4 3; 1 1 4 1 5 2 3 3 2; 3 5 2 2 6 3 1 4 3; 2 6 4 5 1 3 6 5 5; 2 3 6 2 5 6 2 6 3; 4 4 1 6 2 2 6 4 1; 4 3 5 4 1 2 6 1 5; 2 2 1 6 1 4 3 2 1; 3 6 4 2 2 1 4 6 4]]
# starting_configuration_L_11 = [[1 4 3 2 2 4 3 5 2 1 5; 1 3 1 3 1 2 2 1 4 3 1; 1 3 4 5 5 1 3 5 5 2 5; 1 5 6 6 4 2 1 2 6 4 3; 2 2 1 5 2 5 6 6 2 6 5; 2 4 6 5 2 5 4 5 2 3 6; 4 1 6 3 6 6 6 4 5 5 5; 5 4 5 2 3 3 2 2 4 1 6; 4 5 1 6 6 1 6 5 4 2 2; 5 2 2 1 4 1 6 5 2 5 2; 6 3 6 6 2 1 3 5 6 2 4],[4 2 2 1 3 2 5 3 5 3 1; 4 4 6 4 2 6 5 6 1 2 1; 1 1 4 2 3 1 6 2 6 2 6; 3 2 3 6 5 5 4 5 5 4 1; 3 1 2 1 1 3 3 6 3 4 5; 6 2 5 1 6 6 6 3 1 6 3; 5 3 6 5 6 2 1 3 5 3 6; 5 1 4 6 2 5 1 1 1 5 3; 5 2 4 6 4 2 4 2 5 2 5; 4 1 5 4 6 1 4 6 2 2 1; 3 5 2 6 6 1 1 4 5 2 3],[5 2 1 1 1 3 4 5 2 4 3; 6 2 3 4 5 5 2 3 6 6 6; 3 6 5 5 2 1 5 2 4 1 2; 1 2 1 2 3 1 6 4 5 2 5; 6 1 3 4 4 1 2 3 3 4 5; 1 3 5 3 6 3 1 5 1 3 1; 3 6 3 2 1 3 6 5 2 6 1; 4 3 2 3 5 4 2 1 2 5 2; 5 6 1 6 2 6 6 4 5 2 1; 3 1 6 5 5 4 3 5 1 5 4; 6 3 3 1 5 1 4 3 6 5 4],[6 4 1 5 6 3 6 3 2 3 3; 3 5 3 3 6 1 1 6 4 1 1; 1 2 2 5 4 5 3 4 4 1 3; 6 1 6 5 4 5 1 6 6 6 6; 4 3 1 2 4 1 1 5 5 3 6; 4 5 3 3 6 4 6 1 4 4 4; 2 4 6 2 2 1 4 6 3 6 6; 5 4 4 6 1 1 3 4 1 5 4; 4 5 4 4 4 4 3 5 1 6 5; 5 2 4 2 2 5 3 3 3 6 3; 3 2 5 5 5 5 6 5 3 2 6],[6 4 2 2 2 1 3 3 1 1 3; 1 6 4 5 3 3 2 3 6 1 2; 6 4 2 4 3 5 4 2 2 6 2; 3 1 6 5 6 4 6 5 4 2 5; 6 2 4 2 2 3 5 4 1 2 2; 4 6 1 3 3 5 3 3 4 1 2; 3 2 4 1 6 3 4 6 4 6 3; 4 3 5 5 1 1 3 6 1 5 3; 3 3 5 5 6 4 2 1 6 4 1; 1 4 1 6 5 6 3 3 2 1 5; 5 3 2 4 5 2 1 5 2 5 1],[2 3 2 2 1 4 1 2 3 4 6; 6 2 4 4 3 4 1 4 6 4 5; 1 2 3 1 1 1 6 3 1 4 3; 3 1 6 4 6 1 4 3 6 4 3; 3 1 3 4 2 2 6 2 2 1 2; 2 4 6 2 5 3 5 6 3 4 3; 6 2 4 5 4 2 6 1 4 6 4; 4 3 6 5 5 1 3 4 5 4 4; 4 3 6 4 5 2 5 4 5 3 2; 3 1 4 4 5 5 1 1 5 2 6; 2 4 1 4 6 2 6 2 6 4 3]]
# starting_configurations = [starting_configuration_L_5, starting_configuration_L_7, starting_configuration_L_9, starting_configuration_L_11]


experiment_name = "L_9_T_$(sample_temperature)_t_$(autocorrelation_window_length)_trial_$(trial_number)_slice"

annealing_swap_move_probability = 1.0
autocorrelation_swap_move_probability = 0.0 # DOING SLICE-ROTATIONS NOW

fast_autocorrelation_experiment(experiment_name, 9, annealing_swap_move_probability, autocorrelation_swap_move_probability, 10.0, 0.1, 100, sample_temperature, 10000, autocorrelation_window_length; inherent_disorder=true, lag_limit=lag_limit)