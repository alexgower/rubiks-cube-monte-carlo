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





@inbounds @fastmath function neighbour_energy_deltas_distribution_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, neighbour_energy_deltas_sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0, collecting_swap_move_neighbours::Bool=false, neighbour_per_configuration_sample_size::Int64=0, collect_minimum_neighbour_energy_delta_only::Bool=false, extra_swap_moves::Int64=0, extra_slice_rotations::Int64=0)

    neighbour_energy_deltas_sample_temperatures = sort!(neighbour_energy_deltas_sample_temperatures)

    T_1 = 10.0
    T_0 = 0.09
    T_swap = 1.5
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; neighbour_energy_deltas_sample_temperatures]
    sort!(temperature_vector, rev=true)

    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    if relaxation_iterations == 0
        temperature_vector, E_average_by_temperature, _, _, _, _, neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_energy_deltas_sample_temperatures=neighbour_energy_deltas_sample_temperatures, average_sample_size=average_sample_size, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbour_per_configuration_sample_size=neighbour_per_configuration_sample_size, collect_minimum_neighbour_energy_delta_only=collect_minimum_neighbour_energy_delta_only, extra_swap_moves=extra_swap_moves, extra_slice_rotations=extra_slice_rotations)
    else
        relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
        temperature_vector, E_average_by_temperature, _, _, _, _, neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_energy_deltas_sample_temperatures=neighbour_energy_deltas_sample_temperatures, average_sample_size=average_sample_size, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbour_per_configuration_sample_size=neighbour_per_configuration_sample_size, collect_minimum_neighbour_energy_delta_only=collect_minimum_neighbour_energy_delta_only, extra_swap_moves=extra_swap_moves, extra_slice_rotations=extra_slice_rotations)
    end

    # Analyse Results ---------

    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)

        simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)_T=$(sample_temperature)"
        
        # Create plot ----------

        # Find maximum and minima neighbour energy deltas to use as end points of histogram
        max_neighbour_energy_delta = maximum(neighbour_energy_deltas_by_temperature)
        min_neighbour_energy_delta = minimum(neighbour_energy_deltas_by_temperature)

        # Find maximimum probability density
        max_probability_density = 0
        for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
            histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])
            probability_densities = collect(values(histogram))
            probability_densities = probability_densities ./ sum(probability_densities)
            this_max_probability_density = max(max_probability_density, maximum(probability_densities))

            if this_max_probability_density > max_probability_density
                max_probability_density = this_max_probability_density
            end
        end

        try
            
            # Create a histogram
            histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])

            # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
            for energy_delta in min_neighbour_energy_delta:max_neighbour_energy_delta
                if !haskey(histogram, energy_delta)
                    histogram[energy_delta] = 0
                end
            end

            energy_deltas = collect(keys(histogram))
            println("Energy deltas: ", energy_deltas)

            probability_densities = collect(values(histogram))
            println("Probability densities: ", probability_densities)

            # Normalize the values
            probability_densities = probability_densities ./ sum(probability_densities)

            # Plot histogram
            title = collect_minimum_neighbour_energy_delta_only ? "Min Neighbour Energy Deltas Distribution at T=$(sample_temperature)" : "Neighbour Energy Deltas Distribution at T=$(sample_temperature)"
            label = collect_minimum_neighbour_energy_delta_only ? "Min Neighbour Energy Deltas" : "Neighbouring Energy Deltas"
            graph = bar(energy_deltas, probability_densities, xticks=energy_deltas, xtickfontsize=4, label=label, ylim=(0,max_probability_density*1.1), legend=:topleft, title=title, xlabel=label, ylabel="Probability Density")

            # Plot One in a Hundred Energy Dashed Line
            one_in_hundred_energy = log(100)*neighbour_energy_deltas_sample_temperatures[index]
            vline!(graph, [one_in_hundred_energy], linestyle=:dash, label="1/100 Energy", color=:orange)

            # Plot One in a Thousand Energy Dashed Line
            one_in_thousand_energy = log(1000)*neighbour_energy_deltas_sample_temperatures[index]
            vline!(graph, [one_in_thousand_energy], linestyle=:dash, label="1/1000 Energy", color=:red)

            # Save graphs ----------
            savefig(graph, "results/neighbour_energy_deltas_distribution_results/$(simulation_name_to_use).png")


        catch ex

            println("Cannot display or save plot")
            showerror(stdout, ex)

        end

        # Save Neighbour Energy Delta Results ----------
        try
            data_simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)"

            touch(joinpath("results/neighbour_energy_deltas_distribution_results",data_simulation_name_to_use))

            open(joinpath("results/neighbour_energy_deltas_distribution_results",data_simulation_name_to_use), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Trials=$(average_sample_size), Neighbours Per Configuration Sample Size=$(neighbours_per_configuration_sample_size), Collecting Swap Move Neighbours = $(collecting_swap_move_neighbours) \n")
                write(simulation_file, "Sample Temperature T, Average Energy E, Neighbour Energy Deltas \n")
                
                for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
                    write(simulation_file, "$(neighbour_energy_deltas_sample_temperatures[index]),$(E_average_by_temperature[index]), $(neighbour_energy_deltas_by_temperature[index,:]) \n")
                end
            end

        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end

    end


end











@everywhere @inbounds @fastmath function parallel_neighbour_energy_deltas_distribution_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, neighbour_energy_deltas_sample_temperatures::Vector{Float64}, average_sample_size::Int64; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0, collecting_swap_move_neighbours::Bool=false, neighbours_per_configuration_sample_size::Int64=0, collect_minimum_neighbour_energy_delta_only::Bool=false, extra_swap_moves::Int64=0, extra_slice_rotations::Int64=0)

    neighbour_energy_deltas_sample_temperatures = sort!(neighbour_energy_deltas_sample_temperatures)

    T_1 = 10.0
    T_0 = 0.09
    T_swap = 1.5
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; neighbour_energy_deltas_sample_temperatures]
    sort!(temperature_vector, rev=true)

    if neighbours_per_configuration_sample_size == 0
        neighbours_per_configuration_sample_size = configuration_network_degree(L, false)
    end

    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    E_values_by_temperature = SharedArray{Float64,2}(length(temperature_vector), average_sample_size)

    if !collect_minimum_neighbour_energy_delta_only
        neighbour_energy_deltas_by_temperature = SharedArray{Float64,2}(length(neighbour_energy_deltas_sample_temperatures), neighbours_per_configuration_sample_size*average_sample_size)
    else
        neighbour_energy_deltas_by_temperature = SharedArray{Float64,2}(length(neighbour_energy_deltas_sample_temperatures), 1*average_sample_size)
    end

    # neighbour_energy_deltas_by_temperature .= ones(length(neighbour_energy_deltas_sample_temperatures), neighbours_per_configuration_sample_size*average_sample_size)
    # neighbour_energy_deltas_by_temperature = zeros(length(neighbour_energy_deltas_sample_temperatures), neighbours_per_configuration_sample_size*average_sample_size)


    @sync @distributed for trial in 1:average_sample_size
        printstyled("Trial: $trial \n", color=:light_blue)
        if relaxation_iterations == 0
            temperature_vector, trial_E_average_by_temperature, _, _, _, _, trial_neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=false, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_energy_deltas_sample_temperatures=neighbour_energy_deltas_sample_temperatures, average_sample_size=1, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_minimum_neighbour_energy_delta_only=collect_minimum_neighbour_energy_delta_only, extra_swap_moves=extra_swap_moves, extra_slice_rotations=extra_slice_rotations)
        else
            relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
            temperature_vector, trial_E_average_by_temperature, _, _, _, _, trial_neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=false, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_energy_deltas_sample_temperatures=neighbour_energy_deltas_sample_temperatures, average_sample_size=1, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_minimum_neighbour_energy_delta_only=collect_minimum_neighbour_energy_delta_only, extra_swap_moves=extra_swap_moves, extra_slice_rotations=extra_slice_rotations)
        end

        if !collect_minimum_neighbour_energy_delta_only
            neighbour_energy_deltas_by_temperature[:,(trial-1)*neighbours_per_configuration_sample_size+1:trial*neighbours_per_configuration_sample_size] .= trial_neighbour_energy_deltas_by_temperature[:,:]
        else
            neighbour_energy_deltas_by_temperature[:,trial] .= trial_neighbour_energy_deltas_by_temperature[:,:]
        end
        E_values_by_temperature[:,trial] .= trial_E_average_by_temperature[:]
    end

    E_average_by_temperature = [mean(E_values_by_temperature[index,:]) for index in 1:length(temperature_vector)]


    # Analyse Results ---------

    # Find maximum and minima neighbour energy deltas to use as end points of histogram
    max_neighbour_energy_delta = maximum(neighbour_energy_deltas_by_temperature)
    min_neighbour_energy_delta = minimum(neighbour_energy_deltas_by_temperature)
    
    # Find maximimum probability density
    max_probability_density = 0
    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
        histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])
        probability_densities = collect(values(histogram))
        probability_densities = probability_densities ./ sum(probability_densities)
        this_max_probability_density = max(max_probability_density, maximum(probability_densities))

        if this_max_probability_density > max_probability_density
            max_probability_density = this_max_probability_density
        end
    end

    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)

        simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)_T=$(sample_temperature)"
        
        # Create plot ----------

        # Find maximum and minima neighbour energy deltas to use as end points of histogram
        max_neighbour_energy_delta = maximum(neighbour_energy_deltas_by_temperature)
        min_neighbour_energy_delta = minimum(neighbour_energy_deltas_by_temperature)

        # Find maximimum probability density
        max_probability_density = 0
        for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
            histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])
            probability_densities = collect(values(histogram))
            probability_densities = probability_densities ./ sum(probability_densities)
            this_max_probability_density = max(max_probability_density, maximum(probability_densities))

            if this_max_probability_density > max_probability_density
                max_probability_density = this_max_probability_density
            end
        end

        try
            
            # Create a histogram
            histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])

            # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
            for energy_delta in min_neighbour_energy_delta:max_neighbour_energy_delta
                if !haskey(histogram, energy_delta)
                    histogram[energy_delta] = 0
                end
            end

            energy_deltas = collect(keys(histogram))
            println("Energy deltas: ", energy_deltas)

            probability_densities = collect(values(histogram))
            println("Probability densities: ", probability_densities)

            # Normalize the values
            probability_densities = probability_densities ./ sum(probability_densities)

            # Plot histogram
            title = collect_minimum_neighbour_energy_delta_only ? "Min Neighbour Energy Deltas Distribution at T=$(sample_temperature)" : "Neighbour Energy Deltas Distribution at T=$(sample_temperature)"
            label = collect_minimum_neighbour_energy_delta_only ? "Min Neighbour Energy Deltas" : "Neighbouring Energy Deltas"
            graph = bar(energy_deltas, probability_densities, xticks=energy_deltas, xtickfontsize=4, label=label, ylim=(0,max_probability_density*1.1), legend=:topleft, title=title, xlabel=label, ylabel="Probability Density")

            # Plot One in a Hundred Energy Dashed Line
            one_in_hundred_energy = log(100)*neighbour_energy_deltas_sample_temperatures[index]
            vline!(graph, [one_in_hundred_energy], linestyle=:dash, label="1/100 Energy", color=:orange)

            # Plot One in a Thousand Energy Dashed Line
            one_in_thousand_energy = log(1000)*neighbour_energy_deltas_sample_temperatures[index]
            vline!(graph, [one_in_thousand_energy], linestyle=:dash, label="1/1000 Energy", color=:red)

            # Save graphs ----------
            savefig(graph, "results/neighbour_energy_deltas_distribution_results/$(simulation_name_to_use).png")


        catch ex

            println("Cannot display or save plot")
            showerror(stdout, ex)

        end
    end

    # Save Neighbour Energy Delta Results ----------
    try
        data_simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)"

        touch(joinpath("results/neighbour_energy_deltas_distribution_results",data_simulation_name_to_use))

        open(joinpath("results/neighbour_energy_deltas_distribution_results",data_simulation_name_to_use), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Trials=$(average_sample_size), Neighbours Per Configuration Sample Size=$(neighbours_per_configuration_sample_size), Collecting Swap Move Neighbours = $(collecting_swap_move_neighbours) \n")
            write(simulation_file, "Sample Temperature T, Average Energy E, Neighbour Energy Deltas \n")
            
            for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
                write(simulation_file, "$(neighbour_energy_deltas_sample_temperatures[index]),$(E_average_by_temperature[index]), $(neighbour_energy_deltas_by_temperature[index,:]) \n")
            end
        end

    catch ex

        println("Cannot save results to file")
        showerror(stdout, ex)

    end



end






