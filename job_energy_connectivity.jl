number_of_processors = 50
L_values = [9,10]

number_of_L_values = length(L_values)
number_of_processors_per_L_value = Int(floor(number_of_processors/number_of_L_values))


using Distributed
addprocs(number_of_processors)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/neighbour_initial_and_final_energies_distribution_experiment.jl")


@sync @distributed for L_value in L_values
    neighbour_initial_and_final_energies_distribution_experiment("L=$(L_value)_inherent_disorder_E0_E1_slice", L_value, 1.0, 100, [10.0,5.0,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.1]; relaxation_iterations=10000, collecting_swap_move_neighbours=false, neighbours_per_configuration_sample_size=0, average_sample_size_per_temperature=1000, inherent_disorder=true, neighbour_moves_away=1, parallel_anneals=number_of_processors_per_L_value)
end



