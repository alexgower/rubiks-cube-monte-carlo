# Xeon Only ---
# DEPOT_PATH[1]="/u/rscratch/apg59/.julia"
# using Pkg
# Pkg.instantiate()

using Plots
using DelimitedFiles
include("../probes/relaxed_anneal.jl")


# Parallel Case
using Distributed
addprocs(Sys.CPU_THREADS-1)
using SharedArrays
@everywhere include("../probes/relaxed_anneal.jl")







@inbounds @fastmath function energy_histogram_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, energy_histogram_sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0)

    energy_histogram_sample_temperatures = sort!(energy_histogram_sample_temperatures)

    energy_histogram_sample_size = 5000

    T_1 = 10.0
    T_0 = 0.09
    T_swap = 1.5
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; energy_histogram_sample_temperatures]
    sort!(temperature_vector, rev=true)

    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    # TODO remove if want just emergent disorder
    facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
    shuffle!(facelets)
    new_faces = reshape(facelets, 6, L, L)
    for i in 1:6
        cube.configuration[i][:,:] .= new_faces[i,:,:]
    end


    if relaxation_iterations == 0
        temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, measured_relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature, energy_samples_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, energy_histogram_sample_temperatures=energy_histogram_sample_temperatures, average_sample_size=energy_histogram_sample_size)
    else
        relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
        temperature_vector, E_average_by_temperature, E_squared_average_by_temperature, measured_relaxation_iterations_by_temperature, accepted_candidates_by_temperature, final_configuration_correlation_function_by_temperature, energy_samples_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, energy_histogram_sample_temperatures=energy_histogram_sample_temperatures, average_sample_size=energy_histogram_sample_size)
    end

    # Analyse Results ---------

    for (index,sample_temperature) in pairs(energy_histogram_sample_temperatures)

        simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)_T=$(sample_temperature)"
        
        # Create plot ----------

        # Find maximum and minima neighbour energy deltas to use as end points of histogram
        max_energy_sample = maximum(energy_samples_by_temperature)
        min_energy_sample = minimum(energy_samples_by_temperature)

        # Find maximimum probability density
        max_probability_density = 0
        for (index,sample_temperature) in pairs(energy_sample_temperatures)
            histogram = countmap(energy_samples_by_temperature[index,:])
            probability_densities = collect(values(histogram))
            probability_densities = probability_densities ./ sum(probability_densities)
            this_max_probability_density = max(max_probability_density, maximum(probability_densities))

            if this_max_probability_density > max_probability_density
                max_probability_density = this_max_probability_density
            end
        end

        try
            
                
            # Create a histogram
            histogram = countmap(energy_samples_by_temperature[index,:])

            # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
            for energy_sample_value in min_energy_sample:max_energy_sample
                if !haskey(histogram, energy_sample_value)
                    histogram[energy_sample_value] = 0
                end
            end

            energy_sample_values = collect(keys(histogram))
            println("Energy Sample Values: ", energy_sample_values)

            probability_densities = collect(values(histogram))
            println("Probability densities: ", probability_densities)

            # Normalize the values
            probability_densities = probability_densities ./ sum(probability_densities)

            # Plot histogram
            graph = bar(energy_sample_values, probability_densities, xticks=energy_sample_values, xtickfontsize=4, label="Energy", ylim=(0,max_probability_density), legend=:topleft, xlabel="Energy", ylabel="Probability Density")
            L = simulation_name_to_use[3]
            swap_move_probability = simulation_name_to_use[end-3:end]
            title!(graph, "Energy Histogram for L=$L, P_s=$swap_move_probability, T=$sample_temperature")

            # Save graphs ----------
            savefig(graph, "results/energy_histogram_results/$(simulation_name_to_use).png")

        catch ex

            println("Cannot display or save plot")
            showerror(stdout, ex)

        end

        # Save Energy Histogram Results ----------
        try

            touch(joinpath("results/energy_histogram_results",simulation_name_to_use))

            open(joinpath("results/energy_histogram_results",simulation_name_to_use), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, N_sample=$(energy_histogram_sample_size) \n")
                write(simulation_file, "Sample Temperature T, Energy Value \n")
                
                for (index,sample_temperature) in pairs(energy_histogram_sample_temperatures)
                    write(simulation_file, "$(energy_histogram_sample_temperatures[index]), $(energy_samples_by_temperature[index,:]) \n")
                end
            end

        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end

        # Save Relaxed Anneal Results ----------
        try
            normalised_E_average_by_temperature = E_average_by_temperature ./ solved_configuration_energy(cube)
            
            touch(joinpath("results//energy_histogram_results",simulation_name))

            open(joinpath("results/energy_histogram_results",simulation_name), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T \n")
                write(simulation_file, "Temperature T, <E>(T), <-E/E_0>(T), <E^2>(T), Relaxation Iterations=tau(T), Accepted Candidates A(T), Final Configuration Correlation Function Value \n")
                
                for temperature_index in eachindex(E_average_by_temperature)
                    write(simulation_file, "$(temperature_vector[temperature_index]), $(E_average_by_temperature[temperature_index]), $(normalised_E_average_by_temperature[temperature_index]), $(E_squared_average_by_temperature[temperature_index]), $(measured_relaxation_iterations_by_temperature[temperature_index]), $(accepted_candidates_by_temperature[temperature_index]), $(final_configuration_correlation_function_by_temperature[temperature_index]) \n")
                end
            end

        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end
    end


end










@everywhere @inbounds @fastmath function parallel_energy_histogram_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, energy_histogram_sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0)

    energy_histogram_sample_temperatures = sort!(energy_histogram_sample_temperatures)

    # Here just use many trials to exploit parallel speedup (rather than statistical independence)
    # So strike compromise between overhead of many processes and waiting 1*tau in between energy samples in serial case
    sample_size_per_trial = 1500
    number_of_trials = 3
    energy_histogram_sample_size = sample_size_per_trial*number_of_trials

    T_1 = 10.0
    T_0 = 0.09
    T_swap = 1.5
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; energy_histogram_sample_temperatures]
    sort!(temperature_vector, rev=true)

    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    energy_samples_by_temperature = SharedArray{Float64,2}(length(energy_histogram_sample_temperatures), sample_size_per_trial*number_of_trials)

    @sync @distributed for trial in 1:number_of_trials
        printstyled("Trial: $trial \n", color=:light_blue)


        if relaxation_iterations == 0
            temperature_vector, _, _, _, _, _, trial_energy_samples_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, energy_histogram_sample_temperatures=energy_histogram_sample_temperatures, average_sample_size=sample_size_per_trial)
        else
            relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
            temperature_vector, _, _, _, _, _, trial_energy_samples_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, energy_histogram_sample_temperatures=energy_histogram_sample_temperatures, average_sample_size=sample_size_per_trial)
        end

        energy_samples_by_temperature[:,(trial-1)*sample_size_per_trial+1:trial*sample_size_per_trial] .= trial_energy_samples_by_temperature[:,:]

    end

    # Analyse Results ---------

    for (index,sample_temperature) in pairs(energy_histogram_sample_temperatures)

        simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)_T=$(sample_temperature)"
        
        # Create plot ----------

        # Find maximum and minima neighbour energy deltas to use as end points of histogram
        max_energy_sample = maximum(energy_samples_by_temperature)
        min_energy_sample = minimum(energy_samples_by_temperature)

        # Find maximimum probability density
        max_probability_density = 0
        for (index,sample_temperature) in pairs(energy_sample_temperatures)
            histogram = countmap(energy_samples_by_temperature[index,:])
            probability_densities = collect(values(histogram))
            probability_densities = probability_densities ./ sum(probability_densities)
            this_max_probability_density = max(max_probability_density, maximum(probability_densities))

            if this_max_probability_density > max_probability_density
                max_probability_density = this_max_probability_density
            end
        end

        try
            
                
            # Create a histogram
            histogram = countmap(energy_samples_by_temperature[index,:])

            # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
            for energy_sample_value in min_energy_sample:max_energy_sample
                if !haskey(histogram, energy_sample_value)
                    histogram[energy_sample_value] = 0
                end
            end

            energy_sample_values = collect(keys(histogram))
            println("Energy Sample Values: ", energy_sample_values)

            probability_densities = collect(values(histogram))
            println("Probability densities: ", probability_densities)

            # Normalize the values
            probability_densities = probability_densities ./ sum(probability_densities)

            # Plot histogram
            graph = bar(energy_sample_values, probability_densities, xticks=energy_sample_values, xtickfontsize=4, label="Energy", ylim=(0,max_probability_density), legend=:topleft, xlabel="Energy", ylabel="Probability Density")
            L = simulation_name_to_use[3]
            swap_move_probability = simulation_name_to_use[end-3:end]
            title!(graph, "Energy Histogram for L=$L, P_s=$swap_move_probability, T=$sample_temperature")

            # Save graphs ----------
            savefig(graph, "results/energy_histogram_results/$(simulation_name_to_use).png")

        catch ex

            println("Cannot display or save plot")
            showerror(stdout, ex)

        end

        # Save Energy Histogram Results ----------
        try

            touch(joinpath("results/energy_histogram_results",simulation_name_to_use))

            open(joinpath("results/energy_histogram_results",simulation_name_to_use), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, N_sample=$(energy_histogram_sample_size), N_trials=$(number_of_trials), Sample_size_per_trial=$(sample_size_per_trial)  \n")
                write(simulation_file, "Sample Temperature T, Energy Value \n")
                
                for (index,sample_temperature) in pairs(energy_histogram_sample_temperatures)
                    write(simulation_file, "$(energy_histogram_sample_temperatures[index]), $(energy_samples_by_temperature[index,:]) \n")
                end
            end

        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end

    end


end