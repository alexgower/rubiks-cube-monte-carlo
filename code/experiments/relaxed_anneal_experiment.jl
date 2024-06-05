using Plots
using DelimitedFiles

include("../probes/relaxed_anneal.jl")
include("../tools/relaxed_anneal_graphs_plotter.jl")




@inbounds @fastmath function relaxed_anneal_experiment(simulation_name::String, L::Int64, swap_move_probabilities::Vector{Float64}, T_swap::Float64, T_1::Float64, T_0::Float64, N_T::Int64; verbose_metropolis_swap::Bool=false, relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0, bonus_temperatures=[], inherent_disorder::Bool=false, initial_cube_configuration=nothing)

    # Set up temperature vector
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; bonus_temperatures]
    temperature_vector = unique(sort(temperature_vector, rev=true))
    N_T = length(temperature_vector)-1
    
    # Set up arrays to store results
    array_normalised_E_average_by_temperature::Vector{Vector{Float64}} = []
    array_normalised_standard_deviations_by_temperature::Vector{Vector{Float64}} = []
    array_specific_heat_capacities_by_temperature::Vector{Vector{Float64}} = []
    array_M_average_by_temperature::Vector{Vector{Float64}} = []
    array_M_2_average_by_temperature::Vector{Vector{Float64}} = []
    array_M_4_average_by_temperature::Vector{Vector{Float64}} = []


    # Set up cube initial configuration if not provided
    if isnothing(initial_cube_configuration)
        initial_cube_configuration = RubiksCube(L).configuration
        if inherent_disorder
            facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
            shuffle!(facelets)
            new_faces = reshape(facelets, 6, L, L)
            for i in 1:6
                initial_cube_configuration[i][:,:] .= new_faces[i,:,:]
            end
            cube_type = "inherent_disorder_"
        else
            cube_type = "clean_"
        end
    else
        cube_type = "custom_"
    end

    println("Initial Configuration:", initial_cube_configuration)


    # Run experiment for each swap move probability on cube
    for (index,swap_move_probability) in pairs(swap_move_probabilities)
        # Create file name

        simulation_name_to_use = cube_type * simulation_name * '_' * string(swap_move_probability)

        # Create a Rubik's cube object and run annealing function on it
        cube = RubiksCube(L)
        cube.configuration = deepcopy(initial_cube_configuration)

        relaxation_iterations_vector = relaxation_iterations==0 ? nothing : [relaxation_iterations for T in temperature_vector]

        # Run anneal
        temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, M_average_by_temperature, M_2_average_by_temperature, M_4_average_by_temperature, relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap)

        println("Final Configuration:")
        println(cube.configuration)

        # Store normalised results
        push!(array_normalised_E_average_by_temperature, -E_average_by_temperature ./ solved_configuration_energy(cube))
        push!(array_normalised_standard_deviations_by_temperature, sqrt.(E_squared_average_by_temperature .- E_average_by_temperature.^2) ./ solved_configuration_energy(cube))
        push!(array_specific_heat_capacities_by_temperature, (E_squared_average_by_temperature .- E_average_by_temperature.^2) ./ temperature_vector.^2)
        push!(array_M_average_by_temperature, M_average_by_temperature)
        push!(array_M_2_average_by_temperature, M_2_average_by_temperature)
        push!(array_M_4_average_by_temperature, M_4_average_by_temperature)



        # Save Results ----------
        try

            touch(joinpath("results/relaxed_anneal_results",simulation_name_to_use))

            open(joinpath("results/relaxed_anneal_results",simulation_name_to_use), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T \n")
                write(simulation_file, "Temperature T, <E>(T), <-E/E_0>(T), <E^2>(T), c(T), <M>, <M^2>, <M^4> Relaxation Iterations=tau(T), Accepted Candidates A(T), Final Configuration Correlation Function Value \n")
                write(simulation_file, "Initial Configuration: $(initial_cube_configuration) \n")

                for temperature_index in 1:N_T
                    write(simulation_file, "$(temperature_vector[temperature_index]), $(E_average_by_temperature[temperature_index]), $(array_normalised_E_average_by_temperature[index][temperature_index]), $(E_squared_average_by_temperature[temperature_index]), $(array_specific_heat_capacities_by_temperature[index][temperature_index]), $(array_M_average_by_temperature[index][temperature_index]), $(array_M_2_average_by_temperature[index][temperature_index]), $(array_M_4_average_by_temperature[index][temperature_index]), $(relaxation_iterations_by_temperature[temperature_index]), $(accepted_candidates_by_temperature[temperature_index]), $(final_configuration_correlation_function_by_temperature[temperature_index]) \n")
                end

                # write(simulation_file, "\n")
                # write(simulation_file, "# Final Configuration: \n")
                # write(simulation_file, "# $(cube.configuration) \n")
            end

        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end
    end
    

    # PLOT GRAPH
    # relaxed_anneal_graphs_plotter(simulation_name, swap_move_probabilities; inherent_disorder=inherent_disorder)

end


