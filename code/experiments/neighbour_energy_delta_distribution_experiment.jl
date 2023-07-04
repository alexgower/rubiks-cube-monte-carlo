# Xeon Only ---
# DEPOT_PATH[1]="/u/rscratch/apg59/.julia"
# using Pkg
# Pkg.instantiate()

using Plots
using DelimitedFiles

include("../probes/relaxed_anneal.jl")





@inbounds @fastmath function neighbour_energy_deltas_distribution_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, neighbour_energy_deltas_sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0)

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
        temperature_vector, E_average_by_temperature, _, _, _, _, neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_energy_deltas_sample_temperatures=sort!(neighbour_energy_deltas_sample_temperatures), average_sample_size=1000)
    else
        relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
        temperature_vector, E_average_by_temperature, _, _, _, _, neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_energy_deltas_sample_temperatures=sort!(neighbour_energy_deltas_sample_temperatures), average_sample_size=1000)
    end

    # Analyse Results ---------

    for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)

        simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)_T=$(sample_temperature)"
        
        # Create plot ----------

        try
            
            # Create a histogram
            histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])

            # Convert the dictionary to arrays
            energy_deltas = collect(keys(histogram))
            println("Energy deltas: ", energy_deltas)

            probability_densities = collect(values(histogram))
            println("Probability densities: ", probability_densities)

            # Normalize the values
            probability_densities = probability_densities ./ sum(probability_densities)

            # Plot histogram
            graph = bar(energy_deltas, probability_densities, xticks=energy_deltas, xtickfontsize=4, label="Neighbouring Energy Deltas")

            # Plot One in Thousand Energy Dashed Line
            one_in_thousand_energy = log(1000)*neighbour_energy_deltas_sample_temperatures[index]
            vline!(graph, [one_in_thousand_energy], linestyle=:dash, label="1/1000 Energy", color=:red)

            # Save graphs ----------
            savefig(graph, "results/neighbour_energy_deltas_distribution_results/$(simulation_name_to_use).png")

        catch ex

            println("Cannot display or save plot")
            showerror(stdout, ex)

        end

        # Save Results ----------
        try

            touch(joinpath("results/neighbour_energy_deltas_distribution_results",simulation_name_to_use))

            open(joinpath("results/neighbour_energy_deltas_distribution_results",simulation_name_to_use), "w") do simulation_file
                write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T \n")
                write(simulation_file, "Sample Temperature T, Neighbour Energy Deltas \n")
                
                for (index,sample_temperature) in pairs(neighbour_energy_deltas_sample_temperatures)
                    write(simulation_file, "$(neighbour_energy_deltas_sample_temperatures[index]), $(neighbour_energy_deltas_by_temperature[index,:]) \n")
                end
            end

        catch ex

            println("Cannot save results to file")
            showerror(stdout, ex)

        end
    end


end