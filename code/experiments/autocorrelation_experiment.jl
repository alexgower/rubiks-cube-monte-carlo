using Plots
using DelimitedFiles
using LaTeXStrings

include("../probes/autocorrelation_anneal.jl")




@inbounds @fastmath function autocorrelation_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, T_1::Float64, T_0::Float64, N_T::Int64, sample_temperatures::Vector{Float64}, relaxation_iterations_per_temperature::Int64, average_sample_size_per_temperature::Int64, autocorrelation_window_length::Int64; verbose_annealing::Bool=false, verbose_metropolis_swap::Bool=false, original_configuration::Vector{Matrix{Int64}}=empty([[]]))

    cube = RubiksCube(L)

    if !isempty(original_configuration)
        cube.configuration = original_configuration
    end


    # Cover everything in try/except clause
    # try

        temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]

   
        sample_temperatures, energy_autocorrelation_time_by_temperature, energy_stretching_exponent_by_temperature, configuration_autocorrelation_time_by_temperature, configuration_stretching_exponent_by_temperature = autocorrelation_anneal!(cube, temperature_vector, relaxation_iterations_per_temperature, sample_temperatures, average_sample_size_per_temperature, autocorrelation_window_length, swap_move_probability=swap_move_probability, verbose_annealing=verbose_annealing, verbose_metropolis_swap=verbose_metropolis_swap)
   
        println("Finished Annealing")
        println("Energy Autocorrelation Time: $(energy_autocorrelation_time_by_temperature)")
        println("Energy Stretching Exponent: $(energy_stretching_exponent_by_temperature)")
        println("Configuration Autocorrelation Time: $(configuration_autocorrelation_time_by_temperature)")
        println("Configuration Stretching Exponent: $(configuration_stretching_exponent_by_temperature)")

        average_energy_autocorrelation_time_by_temperature = [mean(energy_autocorrelation_time_by_temperature[temperature_index,:]) for temperature_index in 1:length(sample_temperatures)]
        average_energy_stretching_exponent_by_temperature = [mean(energy_stretching_exponent_by_temperature[temperature_index,:]) for temperature_index in 1:length(sample_temperatures)]
        average_configuration_autocorrelation_time_by_temperature = [mean(configuration_autocorrelation_time_by_temperature[temperature_index,:]) for temperature_index in 1:length(sample_temperatures)]
        average_configuration_stretching_exponent_by_temperature = [mean(configuration_stretching_exponent_by_temperature[temperature_index,:]) for temperature_index in 1:length(sample_temperatures)]



        # Save Results ----------

        # try

            # Write  Results to File ----------
            energy_results_filename = simulation_name * '_' * string(swap_move_probability) * ".csv"

            touch(joinpath("results/autocorrelation_anneal_results",energy_results_filename))

            open(joinpath("results/autocorrelation_anneal_results",energy_results_filename), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_1=$T_1, T_0=$T_0, N_T=$N_T, autocorrelatoin_sample_size_per_temperature=$(average_sample_size_per_temperature) ,autocorrelation_window_length=$autocorrelation_window_length \n")
                write(simulation_file, "Original Configuration = $original_configuration \n")
                write(simulation_file, " Sample Temperature T, tau_E samples, beta_E samples, tau_C samples, beta_C samples, tau_E, beta_E, tau_C, beta_C \n")
                
                for temperature_index in 1:length(sample_temperatures)
                    write(simulation_file, "$(sample_temperatures[temperature_index]), $(energy_autocorrelation_time_by_temperature[temperature_index,:]), $(energy_stretching_exponent_by_temperature[temperature_index,:]), $(configuration_autocorrelation_time_by_temperature[temperature_index,:]), $(configuration_stretching_exponent_by_temperature[temperature_index,:]), $(average_energy_autocorrelation_time_by_temperature[temperature_index]), $(average_energy_stretching_exponent_by_temperature[temperature_index]), $(average_configuration_autocorrelation_time_by_temperature[temperature_index]), $(average_configuration_stretching_exponent_by_temperature[temperature_index]) \n")
                end

            end

        # catch ex

        #     println("Cannot save results to file")
        #     showerror(stdout, ex)

        # end
        

        
        # try

            # Create plot of average energy autocorrelation time by temperature and configuration autocorrelation time by temperature
            graph = plot(sample_temperatures, average_energy_autocorrelation_time_by_temperature, xlabel="Temperature", ylabel="Autocorrelation Time, "*L"\tau", title="Rubik's Cube Anneal, L=$L", label="Energy Autocorrelation Time, "*L"\tau_E", marker=(:circle,5))
            savefig(graph, "results/autocorrelation_anneal_results/$(simulation_name)_energy_autocorrelation_time.png")
            display(graph)

            # Create plot of average energy stretching exponent by temperature 
            graph = plot(sample_temperatures, average_energy_stretching_exponent_by_temperature, xlabel="Temperature", ylabel="Stretching Exponent, "*L"\beta", title="Rubik's Cube Anneal, L=$L", label="Energy Stretching Exponent, "*L"\beta_E", marker=(:circle,5))
            savefig(graph, "results/autocorrelation_anneal_results/$(simulation_name)_energy_stretching_exponent.png")
            display(graph)

            # Create plot of average configuration autocorrelation time by temperature
            graph = plot(sample_temperatures, average_configuration_autocorrelation_time_by_temperature, xlabel="Temperature", ylabel="Autocorrelation Time, "*L"\tau", title="Rubik's Cube Anneal, L=$L", label="Configuration Autocorrelation Time, "*L"\tau_C", marker=(:circle,5))
            savefig(graph, "results/autocorrelation_anneal_results/$(simulation_name)_configuration_autocorrelation_time.png")
            display(graph)

            # Create plot of average configuration stretching exponent by temperature
            graph = plot(sample_temperatures, average_configuration_stretching_exponent_by_temperature, xlabel="Temperature", ylabel="Stretching Exponent, "*L"\beta", title="Rubik's Cube Anneal, L=$L", label="Configuration Stretching Exponent, "*L"\beta_C", marker=(:circle,5))
            savefig(graph, "results/autocorrelation_anneal_results/$(simulation_name)_configuration_stretching_exponent.png")
            display(graph)

            # Now plot both on same graph
            graph = plot(sample_temperatures, average_energy_autocorrelation_time_by_temperature, xlabel="Temperature", ylabel="Autocorrelation Time, "*L"\tau", title="Rubik's Cube Anneal, L=$L", label="Energy Autocorrelation Time, "*L"\tau_E", marker=(:circle,5))
            plot!(graph, sample_temperatures, average_configuration_autocorrelation_time_by_temperature, label="Configuration Autocorrelation Time, "*L"\tau_C", marker=(:circle,5))
            savefig(graph, "results/autocorrelation_anneal_results/$(simulation_name)_autocorrelation_time.png")
            display(graph)

            # Now plot both on same graph
            graph = plot(sample_temperatures, average_energy_stretching_exponent_by_temperature, xlabel="Temperature", ylabel="Stretching Exponent, "*L"\beta", title="Rubik's Cube Anneal, L=$L", label="Energy Stretching Exponent, "*L"\beta_E", marker=(:circle,5))
            plot!(graph, sample_temperatures, average_configuration_stretching_exponent_by_temperature, label="Configuration Stretching Exponent, "*L"\beta_C", marker=(:circle,5))
            savefig(graph, "results/autocorrelation_anneal_results/$(simulation_name)_stretching_exponent.png")
            display(graph)

        # catch ex
        #     println("Cannot display or save results")
        #     showerror(stdout, ex)

        # end

    # catch ex
    #     println("General Error")
    #     showerror(stdout, ex)

    #     open("error_log.txt", "a") do error_file
    #         write(error_file, "Some Error: \n" * sprint(showerror(stdout, ex)) * "\n")
    #     end

    # end

end


