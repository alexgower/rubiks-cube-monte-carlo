

using Distributed
number_of_processors = parse(Int, ARGS[2])
addprocs(number_of_processors)

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.instantiate()

@everywhere include("code/experiments/connections_experiment.jl")

starting_configuration_L_5 = [[2 3 3 2 6; 5 3 6 4 4; 1 2 5 1 5; 2 2 6 6 2; 2 1 2 3 2],[1 2 5 5 3; 6 6 1 2 3; 2 5 5 4 3; 4 1 6 5 4; 4 6 6 2 4], [3 4 5 4 6; 2 5 4 4 1; 5 3 2 5 4; 4 6 1 2 5; 2 5 5 4 6], [3 3 1 6 3; 2 3 1 6 3; 6 6 4 3 2; 1 1 2 3 2; 6 3 1 5 3], [2 3 4 4 1; 3 3 1 3 6; 3 4 4 5 1; 4 5 5 4 1; 6 1 1 6 6], [4 1 2 3 5; 2 1 6 6 5; 4 1 5 5 1; 2 1 5 1 6; 5 4 4 3 6]]
starting_configuration_L_7 = [[5 6 4 2 6 3 1; 5 5 4 5 6 5 6; 4 2 2 2 4 1 3; 4 1 4 4 1 1 2; 5 3 1 3 3 2 3; 1 1 1 4 6 5 6; 2 5 3 6 6 2 5],[1 6 2 3 1 2 6; 4 5 3 2 4 3 2; 3 6 4 4 6 3 4; 5 4 1 6 6 5 4; 1 5 2 5 5 3 1; 1 1 2 1 6 2 2; 3 4 2 3 2 6 5], [6 5 5 2 5 5 5; 1 5 4 1 2 3 4; 4 6 5 3 5 6 2; 4 3 1 6 5 2 2; 3 2 6 1 4 5 6; 4 4 6 1 4 5 4; 2 2 3 1 6 4 5], [2 3 2 2 3 3 1; 4 6 5 1 6 4 5; 5 4 6 4 4 3 4; 4 3 1 1 3 3 3; 5 2 1 5 4 4 2; 2 6 5 4 5 6 5; 5 5 6 3 3 1 2], [2 2 2 2 1 3 6; 4 5 3 3 5 5 5; 2 3 1 2 6 4 4; 1 6 1 2 1 6 3; 6 3 6 6 5 1 3; 4 3 1 3 5 1 1; 4 6 3 2 4 4 6], [5 6 3 2 2 3 2; 2 4 2 3 4 5 1; 4 1 1 1 6 3 1; 6 2 4 1 2 1 3; 6 6 6 6 6 3 5; 5 4 2 6 1 3 1; 3 5 3 1 4 1 6]]
starting_configuration_L_9 = [[2 4 2 6 4 5 2 6 2; 6 6 3 5 1 3 4 3 4; 5 6 6 3 1 5 1 4 4; 4 1 5 6 6 3 1 4 1; 2 6 4 4 6 2 4 5 5; 5 5 1 2 2 2 1 4 3; 3 4 1 1 5 3 3 3 6; 3 4 2 3 2 4 3 4 2; 6 6 5 1 1 3 5 4 4], [6 6 6 3 2 6 1 2 4; 5 1 5 1 2 6 6 2 4; 4 5 2 3 3 1 5 1 2; 4 4 3 4 6 1 6 5 1; 5 5 4 6 6 3 4 3 6; 5 5 4 6 2 3 4 1 4; 6 6 6 2 4 6 1 6 6; 5 5 1 5 5 1 5 4 1; 4 3 5 5 4 1 4 2 5], [5 2 5 2 6 3 1 2 2; 5 3 3 4 3 5 5 4 1; 3 5 4 2 3 1 6 5 2; 4 5 3 2 5 4 2 2 3; 5 4 5 5 1 3 2 6 1; 6 4 3 6 5 1 3 1 6; 5 4 6 3 1 5 3 2 1; 6 5 1 4 3 1 3 6 3; 2 4 5 6 2 1 1 3 5], [3 2 3 5 4 6 2 6 3; 4 6 3 4 3 2 3 1 5; 4 1 2 1 2 4 2 5 3; 6 6 6 1 1 2 6 6 5; 2 2 6 6 6 2 4 1 3; 1 1 2 1 4 3 2 1 6; 3 2 1 1 4 1 4 3 2; 5 3 2 3 1 3 2 1 3; 5 5 1 6 1 5 5 1 2], [5 3 4 6 4 2 5 5 6; 4 2 5 5 1 1 5 2 1; 6 1 2 4 6 2 1 1 4; 6 3 1 5 1 6 2 5 6; 1 2 2 3 3 2 5 5 4; 3 6 5 3 1 1 5 4 4; 2 4 5 2 5 3 4 4 6; 2 3 3 3 6 4 5 6 5; 3 3 6 2 1 6 3 4 4], [2 6 4 5 6 3 3 4 3; 1 1 4 1 5 2 3 3 2; 3 5 2 2 6 3 1 4 3; 2 6 4 5 1 3 6 5 5; 2 3 6 2 5 6 2 6 3; 4 4 1 6 2 2 6 4 1; 4 3 5 4 1 2 6 1 5; 2 2 1 6 1 4 3 2 1; 3 6 4 2 2 1 4 6 4]]
starting_configuration_L_11 = [[4 2 1 3 3 2 6 6 6 6 5; 3 3 3 5 5 5 2 4 6 1 6; 3 4 3 6 6 4 6 5 1 6 4; 4 1 2 2 1 2 4 1 2 1 1; 1 5 4 2 5 6 5 4 3 6 3; 6 6 5 2 4 4 5 1 3 6 4; 4 5 2 6 3 2 4 6 6 5 1; 4 6 5 6 6 6 6 5 5 6 4; 5 3 6 1 2 5 4 3 6 5 5; 5 2 6 2 4 6 4 1 6 3 2; 4 3 4 3 3 6 2 3 1 5 1],[4 6 2 4 3 6 1 1 2 2 6; 6 5 4 1 6 2 2 6 4 1 4; 4 1 5 6 1 4 2 2 1 6 6; 6 6 4 6 5 2 6 1 2 5 1; 1 3 2 2 4 1 3 2 4 2 6; 2 4 1 5 2 1 1 2 6 2 5; 2 3 2 4 3 4 4 2 5 6 6; 5 1 4 5 1 3 1 3 5 1 2; 6 3 5 5 4 5 3 5 3 6 4; 3 2 5 6 5 1 6 6 1 3 1; 1 1 5 3 5 6 2 3 1 5 1],[4 1 3 6 1 3 5 2 5 3 6; 1 3 5 2 5 1 3 6 1 2 1; 2 6 5 4 2 2 2 2 6 1 1; 5 1 6 3 2 3 5 3 4 1 5; 6 4 4 5 6 1 5 3 5 4 4; 6 4 1 4 3 2 3 5 1 2 3; 6 3 5 5 3 5 1 3 3 4 5; 2 1 5 5 6 1 6 1 5 4 3; 3 6 5 5 2 4 5 5 3 6 4; 4 5 3 2 2 3 1 4 3 5 4; 2 3 3 1 5 3 6 6 2 6 4],[5 1 3 1 1 4 2 1 3 3 2; 2 1 4 2 4 3 4 3 4 5 2; 1 4 2 2 3 2 1 2 4 6 1; 3 4 3 3 6 3 1 6 3 6 6; 5 3 5 2 4 2 3 1 1 5 4; 3 2 3 5 4 3 1 4 2 2 5; 4 2 6 6 4 3 4 6 4 3 1; 6 3 1 5 5 5 6 2 6 4 5; 3 4 3 3 6 1 1 5 2 2 2; 5 2 4 4 3 2 3 3 1 3 3; 2 3 4 2 2 3 1 2 2 2 2],[4 6 4 2 3 4 6 6 5 2 4; 1 2 4 6 5 6 3 3 3 6 3; 5 4 4 5 3 1 1 4 2 6 3; 1 6 5 2 3 6 2 4 2 6 6; 5 3 2 6 6 1 6 2 6 4 5; 5 6 2 6 4 5 6 5 6 5 1; 4 5 3 2 1 1 6 4 1 1 1; 5 6 3 3 2 3 1 3 6 2 3; 2 3 6 1 4 4 5 4 3 6 3; 3 1 2 4 2 5 6 4 1 6 6; 1 2 3 2 4 6 1 1 5 3 6],[1 1 1 1 5 2 6 1 5 3 5; 4 3 2 2 5 5 4 3 3 4 3; 4 4 1 5 2 3 5 5 2 5 5; 2 6 2 6 2 1 6 1 3 1 2; 3 5 3 5 4 2 6 4 4 4 4; 1 1 1 5 2 5 4 3 1 4 2; 4 4 2 4 1 4 1 1 2 6 5; 2 1 1 4 1 5 1 4 6 2 5; 5 5 1 4 5 5 4 1 5 1 5; 4 4 3 1 2 6 4 3 3 2 4; 2 1 5 6 5 1 6 5 3 4 4]]

L_value = parse(Int, ARGS[1])
extra_save_name = ARGS[3]
println("Extra save name: $extra_save_name")

# L_values = [11]
# starting_configurations = [starting_configuration_L_5, starting_configuration_L_7, starting_configuration_L_9, starting_configuration_L_11]
# starting_configurations[Int((L_value-3)/2)]


sample_temperatures = [10.0*(0.1/10.0)^(m/10) for m in 0:10]
special_sample_temperature = collect(LinRange(0.8,1.5,50))
sample_temperatures = vcat(sample_temperatures, special_sample_temperature)
sort!(sample_temperatures, rev=true)

connections_experiment("$(extra_save_name)_disorder_average_connections_L=$(L_value)_inherent_disorder_E0_E1_E2_slice", L_value, 1.0, 5, 
        sample_temperatures; relaxation_iterations=10000, connections_to_measure="slice", 
        connections_per_configuration_sample_size=0, neighbour_order_to_measure_to=2, average_sample_size_per_temperature=2, 
        initial_cube_configuration=nothing, parallel_anneals=number_of_processors, disorder_average=true, cumulative_connections_measurement=true)



