# Xeon Only ---
# DEPOT_PATH[1]="/u/rscratch/apg59/.julia"
# using Pkg
# Pkg.instantiate()

using Plots
using DelimitedFiles

include("../probes/relaxed_anneal.jl")





@inbounds @fastmath function relaxed_anneal_experiment(simulation_name::String, L::Int64, swap_move_probabilities::Vector{Float64}, T_swap::Float64, T_1::Float64, T_0::Float64, N_T::Int64; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0)


    # Cover everything in try/except clause so can print errors to file if running remotely
    try
        temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
        normalised_E_average_by_temperature = []
        infinite_temperature_normalised_E_average_by_temperature = []


        for (index,swap_move_probability) in pairs(swap_move_probabilities)
            simulation_name_to_use = simulation_name * '_' * string(swap_move_probability)



            # Run Rubik's Cube Anneal ----------

            # Create a Rubik's cube object and run annealing function on it
            cube = RubiksCube(L)

            if relaxation_iterations == 0
                temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap)
            else
                relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
                temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap)
            end

            println("Final Configuration:")
            println(cube.configuration)

            push!(normalised_E_average_by_temperature, -E_average_by_temperature ./ solved_configuration_energy(cube))
            push!(infinite_temperature_normalised_E_average_by_temperature, -E_average_by_temperature ./ infinite_temperature_energy(cube))


            # Save Results ----------
            try
                touch(joinpath("results//relaxed_anneal_results",simulation_name_to_use))

                open(joinpath("results/relaxed_anneal_results",simulation_name_to_use), "w") do simulation_file
                    write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T \n")
                    write(simulation_file, "Temperature T, <E>(T), <-E/E_0>(T), <E^2>(T), Relaxation Iterations=tau(T), Accepted Candidates A(T), Final Configuration Correlation Function Value \n")
                    
                    for temperature_index in 1:N_T+1
                        write(simulation_file, "$(temperature_vector[temperature_index]), $(E_average_by_temperature[temperature_index]), $(normalised_E_average_by_temperature[index][temperature_index]), $(E_squared_average_by_temperature[temperature_index]), $(relaxation_iterations_by_temperature[temperature_index]), $(accepted_candidates_by_temperature[temperature_index]), $(final_configuration_correlation_function_by_temperature[temperature_index]) \n")
                    end
                end

            catch ex

                println("Cannot save results to file")
                showerror(stdout, ex)

            end
        end
        

        
        try
            # Create plot ----------
            
            if normalization == "solved"

                graph = plot(temperature_vector, normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Solved Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
                hline!(graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
                hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")

            else # normalization == "infinite_temperature" case

                graph = plot(temperature_vector, infinite_temperature_normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Infinite Temperature Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
                hline!(graph, [-1.0], linestyle=:dash, color=:black, label="")
                hline!(graph, [-6.0], linestyle=:dash, color=:black, label="")

            end

            vline!(graph, [T_swap], linestyle=:dash, label="Swap Moves Enabled")

            # Add other data to graph for comparison if exists ----------
            if isfile(joinpath("results/relaxed_anneal_results","other_data.csv"))

                data_matrix = readdlm(joinpath("results/relaxed_anneal_results","other_data.csv"), ',', Float64, '\n', skipstart=0)

                other_temperature_vector = copy(data_matrix[:,1])
                other_data = data_matrix[:,2]
            
                plot!(graph, other_temperature_vector, other_data, label="Ollie Results", seriestype=:scatter, color="blue", ms=2, ma=0.5)

            end

            # Save graph ----------
            savefig(graph, "results//relaxed_anneal_results/$simulation_name.png")

        catch ex

            println("Cannot display or save results")
            showerror(stdout, ex)

        end

    catch ex

        println("General Error")
        showerror(stdout, ex)

    end
end


