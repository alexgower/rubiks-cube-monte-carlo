#--- To Do: ---

# - possibly integrate energy/order parameter histories in MC program
# - add average order parameter and order parameter squared alongside average energy measurement
# - use ensemble work

# ---

using Plots

include("rubiks_cube.jl")
include("monte_carlo.jl")
include("anneal.jl")
include("swap_moves.jl")

@inbounds @fastmath function swap_prob_anneal_experiment(simulation_name::String, L::Int64, swap_move_probabilities::Vector{Float64}, T_swap::Float64, T_1::Float64, T_0::Float64, N_T::Int64; verbose_metropolis_swap::Bool=false, normalization::String="solved")


    # Cover everything in try/except clause so can print errors to file if running remotely
    try
        temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
        normalised_E_average_by_temperature = []
        infinite_temperature_normalised_E_average_by_temperature = []
        combined_relaxation_iterations_by_temperature = []

        for (index,swap_move_probability) in pairs(swap_move_probabilities)
            simulation_name_to_use = simulation_name * '_' * string(swap_move_probability)

            # Run Rubik's Cube Anneal ----------

            # Create a Rubik's cube object and run annealing function on it
            cube = RubiksCube(L)
            temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature = anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_finder_mode=true)

            println("Final Configuration:")
            println(cube.configuration)

            push!(normalised_E_average_by_temperature, -E_average_by_temperature ./ solved_configuration_energy(cube))
            push!(infinite_temperature_normalised_E_average_by_temperature, -E_average_by_temperature ./ infinite_temperature_energy(cube))
            push!(combined_relaxation_iterations_by_temperature, relaxation_iterations_by_temperature)


            # Save Results ----------
            try
                touch(joinpath("results",simulation_name_to_use))

                open(joinpath("results",simulation_name_to_use), "w") do simulation_file
                    write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T \n")
                    write(simulation_file, "Temperature T, <E>(T), <-E/E_0>(T), <E^2>(T), Relaxation Iterations=tau(T), Accepted Candidates A(T), Final Configuration Correlation Function Value \n")
                    
                    for temperature_index in 1:N_T+1
                        write(simulation_file, "$(temperature_vector[temperature_index]), $(E_average_by_temperature[temperature_index]), $(normalised_E_average_by_temperature[index][temperature_index]), $(E_squared_average_by_temperature[temperature_index]), $(relaxation_iterations_by_temperature[temperature_index]), $(accepted_candidates_by_temperature[temperature_index]), $(final_configuration_correlation_function_by_temperature[temperature_index]) \n")
                    end
                end

            catch ex
                println("Cannot save results to file")
                showerror(stdout, ex)

                open("error_log.txt", "a") do error_file
                    write(error_file, "Cannot save results to file: \n" * sprint(showerror(stdout, ex)) * "\n")
                end
            end
        end
        

        # Display and Save Plot ----------
        try
            if normalization == "solved"
                plot(temperature_vector, normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Solved Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
            elseif normalization == "infinite_temperature"
                plot(temperature_vector, alternative_normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Infinite Temperature Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities),1,length(swap_move_probabilities)))
            end

            vline!([T_swap], linestyle=:dash, label="Swap Moves Enabled")

            savefig("results/$simulation_name.png")

            # Relaxation iterations plotter ---
            plot(temperature_vector, combined_relaxation_iterations_by_temperature, xlabel="Temperature", ylabel="Relaxation Iterations", title="Relaxation Iterations by Temperature", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
            vline!([T_swap], linestyle=:dash, label="Swap Moves Enabled")

        catch ex
            println("Cannot display or save results")
            showerror(stdout, ex)

            open("error_log.txt", "a") do error_file
                write(error_file, "Cannot save results to file: \n" * sprint(showerror(stdout, ex)) * "\n")
            end
        end

    catch ex
        showerror(stdout, ex)

        open("error_log.txt", "a") do error_file
            write(error_file, "Some Error: \n" * sprint(showerror(stdout, ex)) * "\n")
        end

    end
end

