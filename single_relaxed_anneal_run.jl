# Command-line arguments
model = ARGS[1]
L = parse(Int, ARGS[2])
swap_move_probability = parse(Float64, ARGS[3])
trial_number = parse(Int, ARGS[4])

# Print out the command-line arguments
println("Model: ", model)
println("L: ", L)
println("Swap move probability: ", swap_move_probability)
println("Trial number: ", trial_number)


# Set up packages
println("Instantiating packages")
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Include code file
println("Including experiment file")
include("code/experiments/relaxed_anneal_experiment.jl")



println("Running experiment")

# experiment_name = "L_$(L)_trial_$(trial_number)"
experiment_name = "test"

T_swap = 3.0
T_1 = 10.0
T_0 = 0.3
relaxation_iterations = 10000


if model == "clean"
    inherent_disorder = false
    initial_cube_configuration = nothing

    N_T = 100
    bonus_temperatures = collect(LinRange(0.6,1.0,100))
elseif model == "inherent_disorder"
    inherent_disorder = true
    initial_cube_configuration = nothing

    N_T = 80
    bonus_temperatures = collect(LinRange(0.6,1.0,20))

elseif model == "custom"
    inherent_disorder = false

    N_T = 80
    bonus_temperatures = collect(LinRange(0.6,1.0,20))

    starting_configuration_L_5 = [[2 3 3 2 6; 5 3 6 4 4; 1 2 5 1 5; 2 2 6 6 2; 2 1 2 3 2],[1 2 5 5 3; 6 6 1 2 3; 2 5 5 4 3; 4 1 6 5 4; 4 6 6 2 4], [3 4 5 4 6; 2 5 4 4 1; 5 3 2 5 4; 4 6 1 2 5; 2 5 5 4 6], [3 3 1 6 3; 2 3 1 6 3; 6 6 4 3 2; 1 1 2 3 2; 6 3 1 5 3], [2 3 4 4 1; 3 3 1 3 6; 3 4 4 5 1; 4 5 5 4 1; 6 1 1 6 6], [4 1 2 3 5; 2 1 6 6 5; 4 1 5 5 1; 2 1 5 1 6; 5 4 4 3 6]]
    starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 4 1 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 6 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 2; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 6 5 2 2; 3 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 4 3 1 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 2 1 6 3; 6 3 6 6 5 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
    starting_configuration_L_9 = [[2 4 2 6 4 5 2 6 2; 6 6 3 5 1 3 4 3 4; 5 6 6 3 1 5 1 4 4; 4 1 5 6 6 3 1 4 1; 2 6 4 4 6 2 4 5 5; 5 5 1 2 2 2 1 4 3; 3 4 1 1 5 3 3 3 6; 3 4 2 3 2 4 3 4 2; 6 6 5 1 1 3 5 4 4], [6 6 6 3 2 6 1 2 4; 5 1 5 1 2 6 6 2 4; 4 5 2 3 3 1 5 1 2; 4 4 3 4 6 1 6 5 1; 5 5 4 6 6 3 4 3 6; 5 5 4 6 2 3 4 1 4; 6 6 6 2 4 6 1 6 6; 5 5 1 5 5 1 5 4 1; 4 3 5 5 4 1 4 2 5], [5 2 5 2 6 3 1 2 2; 5 3 3 4 3 5 5 4 1; 3 5 4 2 3 1 6 5 2; 4 5 3 2 5 4 2 2 3; 5 4 5 5 1 3 2 6 1; 6 4 3 6 5 1 3 1 6; 5 4 6 3 1 5 3 2 1; 6 5 1 4 3 1 3 6 3; 2 4 5 6 2 1 1 3 5], [3 2 3 5 4 6 2 6 3; 4 6 3 4 3 2 3 1 5; 4 1 2 1 2 4 2 5 3; 6 6 6 1 1 2 6 6 5; 2 2 6 6 6 2 4 1 3; 1 1 2 1 4 3 2 1 6; 3 2 1 1 4 1 4 3 2; 5 3 2 3 1 3 2 1 3; 5 5 1 6 1 5 5 1 2], [5 3 4 6 4 2 5 5 6; 4 2 5 5 1 1 5 2 1; 6 1 2 4 6 2 1 1 4; 6 3 1 5 1 6 2 5 6; 1 2 2 3 3 2 5 5 4; 3 6 5 3 1 1 5 4 4; 2 4 5 2 5 3 4 4 6; 2 3 3 3 6 4 5 6 5; 3 3 6 2 1 6 3 4 4], [2 6 4 5 6 3 3 4 3; 1 1 4 1 5 2 3 3 2; 3 5 2 2 6 3 1 4 3; 2 6 4 5 1 3 6 5 5; 2 3 6 2 5 6 2 6 3; 4 4 1 6 2 2 6 4 1; 4 3 5 4 1 2 6 1 5; 2 2 1 6 1 4 3 2 1; 3 6 4 2 2 1 4 6 4]]
    starting_configuration_L_11 = [[1 4 3 2 2 4 3 5 2 1 5; 1 3 1 3 1 2 2 1 4 3 1; 1 3 4 5 5 1 3 5 5 2 5; 1 5 6 6 4 2 1 2 6 4 3; 2 2 1 5 2 5 6 6 2 6 5; 2 4 6 5 2 5 4 5 2 3 6; 4 1 6 3 6 6 6 4 5 5 5; 5 4 5 2 3 3 2 2 4 1 6; 4 5 1 6 6 1 6 5 4 2 2; 5 2 2 1 4 1 6 5 2 5 2; 6 3 6 6 2 1 3 5 6 2 4],[4 2 2 1 3 2 5 3 5 3 1; 4 4 6 4 2 6 5 6 1 2 1; 1 1 4 2 3 1 6 2 6 2 6; 3 2 3 6 5 5 4 5 5 4 1; 3 1 2 1 1 3 3 6 3 4 5; 6 2 5 1 6 6 6 3 1 6 3; 5 3 6 5 6 2 1 3 5 3 6; 5 1 4 6 2 5 1 1 1 5 3; 5 2 4 6 4 2 4 2 5 2 5; 4 1 5 4 6 1 4 6 2 2 1; 3 5 2 6 6 1 1 4 5 2 3],[5 2 1 1 1 3 4 5 2 4 3; 6 2 3 4 5 5 2 3 6 6 6; 3 6 5 5 2 1 5 2 4 1 2; 1 2 1 2 3 1 6 4 5 2 5; 6 1 3 4 4 1 2 3 3 4 5; 1 3 5 3 6 3 1 5 1 3 1; 3 6 3 2 1 3 6 5 2 6 1; 4 3 2 3 5 4 2 1 2 5 2; 5 6 1 6 2 6 6 4 5 2 1; 3 1 6 5 5 4 3 5 1 5 4; 6 3 3 1 5 1 4 3 6 5 4],[6 4 1 5 6 3 6 3 2 3 3; 3 5 3 3 6 1 1 6 4 1 1; 1 2 2 5 4 5 3 4 4 1 3; 6 1 6 5 4 5 1 6 6 6 6; 4 3 1 2 4 1 1 5 5 3 6; 4 5 3 3 6 4 6 1 4 4 4; 2 4 6 2 2 1 4 6 3 6 6; 5 4 4 6 1 1 3 4 1 5 4; 4 5 4 4 4 4 3 5 1 6 5; 5 2 4 2 2 5 3 3 3 6 3; 3 2 5 5 5 5 6 5 3 2 6],[6 4 2 2 2 1 3 3 1 1 3; 1 6 4 5 3 3 2 3 6 1 2; 6 4 2 4 3 5 4 2 2 6 2; 3 1 6 5 6 4 6 5 4 2 5; 6 2 4 2 2 3 5 4 1 2 2; 4 6 1 3 3 5 3 3 4 1 2; 3 2 4 1 6 3 4 6 4 6 3; 4 3 5 5 1 1 3 6 1 5 3; 3 3 5 5 6 4 2 1 6 4 1; 1 4 1 6 5 6 3 3 2 1 5; 5 3 2 4 5 2 1 5 2 5 1],[2 3 2 2 1 4 1 2 3 4 6; 6 2 4 4 3 4 1 4 6 4 5; 1 2 3 1 1 1 6 3 1 4 3; 3 1 6 4 6 1 4 3 6 4 3; 3 1 3 4 2 2 6 2 2 1 2; 2 4 6 2 5 3 5 6 3 4 3; 6 2 4 5 4 2 6 1 4 6 4; 4 3 6 5 5 1 3 4 5 4 4; 4 3 6 4 5 2 5 4 5 3 2; 3 1 4 4 5 5 1 1 5 2 6; 2 4 1 4 6 2 6 2 6 4 3]]
    starting_configurations = [starting_configuration_L_5, starting_configuration_L_7, starting_configuration_L_9, starting_configuration_L_11]

    initial_cube_configuration = starting_configurations[Int((L-3)/2)]
else
    println("Invalid model")
    exit(1)
end

relaxed_anneal_experiment(experiment_name, L, [swap_move_probability], T_swap, T_1, T_0, N_T; verbose_metropolis_swap=false, relaxation_iterations=relaxation_iterations, mixing_p_swap=0.0, bonus_temperatures=bonus_temperatures, inherent_disorder=inherent_disorder, initial_cube_configuration=initial_cube_configuration)
