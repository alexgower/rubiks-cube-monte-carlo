using DelimitedFiles
using Plots
using LaTeXStrings
using StatsBase
using CSV
using DataFrames


using GMT

include("../core/rubiks_cube.jl")
include("../core/swap_moves.jl")

function relaxed_anneal_graphs_plotter(simulation_name::String, swap_move_probabilities::Vector{Float64}; normalization::String="solved")

    ### IMPORT DATA

    ### --- SET UP DEFAULT PARAMETERS ---
    header_line = readlines(joinpath("results/relaxed_anneal_results/"*simulation_name*"_1.0"))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    temperature_vector::Vector{Float64} = []
    array_normalised_E_average_by_temperature::Vector{Vector{Float64}} = []
    array_normalised_standard_deviations_by_temperature::Vector{Vector{Float64}} = []
    array_specific_heat_capacities_by_temperature::Vector{Vector{Float64}} = []

    for (index,swap_move_probability) in pairs(swap_move_probabilities)
        simulation_name_to_use = simulation_name * '_' * string(swap_move_probability)

        ### --- READ IN THE DATA ---
        data_matrix = readdlm(joinpath("results/relaxed_anneal_results",simulation_name_to_use), ',', Float64, '\n', skipstart=2)

        temperature_vector = copy(data_matrix[:,1])
        E_average_by_temperature = copy(data_matrix[:,2]) 
        E_squared_average_by_temperature = copy(data_matrix[:,4]) 
        normalization_energy = solved_configuration_energy(cube)

        push!(array_normalised_E_average_by_temperature, -E_average_by_temperature ./ normalization_energy)
        push!(array_normalised_standard_deviations_by_temperature, sqrt.(E_squared_average_by_temperature .- E_average_by_temperature.^2) ./ normalization_energy)
        push!(array_specific_heat_capacities_by_temperature, (E_squared_average_by_temperature .- E_average_by_temperature.^2) ./ temperature_vector.^2)

        # ------------------------------
    end

    

    ### PLOT AND SAVE graphs
    N_T = length(temperature_vector)-1
    normalization_string = normalization=="solved" ? "Solved" : "Infinite Temperature"
            
    mean_std_graph = plot(temperature_vector, array_normalised_E_average_by_temperature, yerr=transpose(array_normalised_standard_deviations_by_temperature), markerstrokecolor=:auto, xlabel="Temperature", ylabel="-Average Energy/"*normalization_string*" Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))
    mean_graph = plot(temperature_vector, array_normalised_E_average_by_temperature, xlabel="Temperature", ylabel="-Average Energy/Solved Energy", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))

    hline!(mean_std_graph, [-1.0], linestyle=:dash, color=:black, label="")
    hline!(mean_graph, [-1.0], linestyle=:dash, color=:black, label="")

    if normalization == "solved"
        hline!(mean_std_graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
        hline!(mean_graph, [-0.16666666666666666], linestyle=:dash, color=:black, label="")
    else # normalization == "infinite_temperature" case
        hline!(mean_std_graph, [-6.0], linestyle=:dash, color=:black, label="")
        hline!(mean_graph, [-6.0], linestyle=:dash, color=:black, label="")
    end

    # Add other data to mean graph for comparison if exists ----------
    if isfile(joinpath("results/relaxed_anneal_results","other_data.csv"))

        data_matrix = readdlm(joinpath("results/relaxed_anneal_results","other_data.csv"), ',', Float64, '\n', skipstart=0)

        other_temperature_vector = copy(data_matrix[:,1])
        other_data = data_matrix[:,2]
    
        plot!(graph, other_temperature_vector, other_data, label="Ollie Results", seriestype=:scatter, color="blue", ms=2, ma=0.5)

    end

    # Create specific heat capacity/temperature plot ----------
    specific_heat_capacity_graph = plot(temperature_vector, array_specific_heat_capacities_by_temperature, xlabel="Temperature", ylabel="Specific Heat Capacity", title="Rubik's Cube Anneal, L=$L", labels=reshape(["P_swap = $swap_move_probability" for swap_move_probability in swap_move_probabilities],1,length(swap_move_probabilities)))











    ### MESSY ENTROPY PLOTTING STUFF


    # Create entropy against energy plot --------
    swap_move_one_probabilty_index = findfirst(swap_move_probabilities .== 1.0)
    
    entropy_by_temperature = zeros(N_T)
    alternative_entropy_by_temperature = zeros(N_T)

    temperatures = reverse(temperature_vector)
    heat_capacities_by_temperature = reverse(array_specific_heat_capacities_by_temperature[swap_move_one_probabilty_index])
    normalized_energies_by_temperature = reverse(array_normalised_E_average_by_temperature[swap_move_one_probabilty_index])
    absolute_energies_by_temperature = (normalized_energies_by_temperature .+ 1.0)*-solved_configuration_energy(cube)

    for temperature_index in 2:N_T
        entropy_by_temperature[temperature_index] = entropy_by_temperature[temperature_index-1] + 0.5 * (temperatures[temperature_index] - temperatures[temperature_index-1]) * (heat_capacities_by_temperature[temperature_index]/temperatures[temperature_index] + heat_capacities_by_temperature[temperature_index-1]/temperatures[temperature_index-1])
    end

    for temperature_index in 2:N_T
        alternative_entropy_by_temperature[temperature_index] = alternative_entropy_by_temperature[temperature_index-1] + 0.5 * (absolute_energies_by_temperature[temperature_index] - absolute_energies_by_temperature[temperature_index-1]) * (1/temperatures[temperature_index] + 1/temperatures[temperature_index-1])
    end

    # Plot entropy against temperature
    entropy_by_temperature_graph = plot(temperatures[2:end], entropy_by_temperature, xlabel="Temperature", ylabel="Entropy", title="Rubik's Cube Anneal, L=$L", label="Equilibrium Entropy", color=:black)
    plot!(entropy_by_temperature_graph, temperatures[2:end], alternative_entropy_by_temperature, label="Alternative Entropy", color=:blue)





    # Plot entropy against absolute energy
    entropy_by_absolute_energy_graph = plot(absolute_energies_by_temperature[2:end], entropy_by_temperature, xlabel="Average Energy from Solved Configuration Energy", ylabel="Entropy", title="Rubik's Cube Anneal, L=$L"; label="Equilibrium Entropy", color=:black, legend=:topleft)

    alternative_entropy_by_absolute_energy_graph = plot(absolute_energies_by_temperature[2:end], alternative_entropy_by_temperature, xlabel="Average Energy from Solved Configuration Energy", ylabel="Entropy", title="Rubik's Cube Anneal, L=$L"; label="Entropy", color=:black, legend=:topleft)
    
    
    slice_rotation_entropy_scaling_gradient = log(6(L-1))*(1/(2*L))
    plot!(entropy_by_absolute_energy_graph, [minimum(absolute_energies_by_temperature[2:end]), maximum(absolute_energies_by_temperature[2:end])], [0, slice_rotation_entropy_scaling_gradient*(maximum(absolute_energies_by_temperature[2:end])-minimum(absolute_energies_by_temperature[2:end]))], line=:dash, color=:blue, lw=2, label=L"\Delta E = 2L"*" (Mode) Layers")
    slice_rotation_entropy_scaling_gradient = log(6(L-1))*(1/(8*L))
    plot!(entropy_by_absolute_energy_graph, [minimum(absolute_energies_by_temperature[2:end]), maximum(absolute_energies_by_temperature[2:end])], [0, slice_rotation_entropy_scaling_gradient*(maximum(absolute_energies_by_temperature[2:end])-minimum(absolute_energies_by_temperature[2:end]))], line=:dash, color=:green, lw=2, label=L"\Delta E = 8L"*" (Max) Layers")


    # Reconstruct (E_0, E_1) histogram from data ----------
    normalised_E0_E1_histogram = reconstruct_histogram(joinpath("results/relaxed_anneal_results",simulation_name*"_normalised_histogram_data.csv"), cube)

    saddle_configurations_from_below_by_absolute_energies = zeros(N_T)
    saddle_configurations_from_even_by_absolute_energies = zeros(N_T)

    alternative_saddle_configurations_from_below_by_absolute_energies = zeros(N_T)
    alternative_saddle_configurations_from_even_by_absolute_energies = zeros(N_T)

    for (index,absolute_energy) in pairs(absolute_energies_by_temperature[2:end])
        # First find histogram bin corresponding to absolute energy
        E1_bin_index = findfirst(>(absolute_energy), normalised_E0_E1_histogram.edges[2]) -1


        for E0_bin_index in 1:E1_bin_index-1

            # TODO checks
            z_slice = 6*(L-1) 
            z_swap = 5940848633430244 # for L=7
            # z_swap = 281635556 # for L=5
            z = z_slice

            # Next sum over all E0 with connections from below to this E_1
            proportion_of_E0_connections_to_E1 = normalised_E0_E1_histogram.weights[E0_bin_index, E1_bin_index]

            E0_index_in_absolute_energies_by_temperature = findfirst(>=(normalised_E0_E1_histogram.edges[1][E0_bin_index]), absolute_energies_by_temperature[2:end])
            configurations_of_E0 = exp(entropy_by_temperature[E0_index_in_absolute_energies_by_temperature])

            saddle_configurations_from_below_by_absolute_energies[index] += proportion_of_E0_connections_to_E1 * configurations_of_E0 * z
            alternative_saddle_configurations_from_below_by_absolute_energies[index] += proportion_of_E0_connections_to_E1 * exp(alternative_entropy_by_temperature[E0_index_in_absolute_energies_by_temperature]) * z


            # Next sum from this E_1 to below
            proportion_of_E1_connections_to_E0 = normalised_E0_E1_histogram.weights[E1_bin_index, E0_bin_index]

            E1_index_in_absolute_energies_by_temperature = findfirst(>=(normalised_E0_E1_histogram.edges[2][E1_bin_index]), absolute_energies_by_temperature[2:end])
            configurations_of_E1 = exp(entropy_by_temperature[E1_index_in_absolute_energies_by_temperature])

            saddle_configurations_from_even_by_absolute_energies[index] += proportion_of_E1_connections_to_E0 * configurations_of_E1 * z
            alternative_saddle_configurations_from_even_by_absolute_energies[index] += proportion_of_E1_connections_to_E0 * exp(alternative_entropy_by_temperature[E1_index_in_absolute_energies_by_temperature]) * z

        end

    end

    # Print raw configuration results
    for i in 1:length(saddle_configurations_from_below_by_absolute_energies)-1
        println("Saddle Configurations at $(absolute_energies_by_temperature[i+1]): ", saddle_configurations_from_below_by_absolute_energies[i])
        println("Alternative Saddle Configurations at $(absolute_energies_by_temperature[i+1]): ", alternative_saddle_configurations_from_below_by_absolute_energies[i])
        println("Total Configurations at $(absolute_energies_by_temperature[i+1]): ", exp(entropy_by_temperature[i+1]))
    end


    entropy_of_saddle_configurations_from_below_by_absolute_energies = log.(saddle_configurations_from_below_by_absolute_energies)
    entropy_of_saddle_configurations_from_even_by_absolute_energies = log.(saddle_configurations_from_even_by_absolute_energies)
    
    plot!(entropy_by_absolute_energy_graph, absolute_energies_by_temperature[2:end], entropy_of_saddle_configurations_from_below_by_absolute_energies, label="Saddle Configurations From Below", color=:red, lw=2)
    plot!(entropy_by_absolute_energy_graph, absolute_energies_by_temperature[2:end], entropy_of_saddle_configurations_from_even_by_absolute_energies, label="Saddle Configurations From Even", color=:orange, lw=2)


    alternative_entropy_of_saddle_configurations_from_below_by_absolute_energies = log.(alternative_saddle_configurations_from_below_by_absolute_energies)
    alternative_entropy_of_saddle_configurations_from_even_by_absolute_energies = log.(alternative_saddle_configurations_from_even_by_absolute_energies)

    plot!(alternative_entropy_by_absolute_energy_graph, absolute_energies_by_temperature[2:end], alternative_entropy_of_saddle_configurations_from_below_by_absolute_energies, label="Alternative Saddle Configurations From Below", color=:red, lw=2)
    plot!(alternative_entropy_by_absolute_energy_graph, absolute_energies_by_temperature[2:end], alternative_entropy_of_saddle_configurations_from_even_by_absolute_energies, label="Alternative Saddle Configurations From Even", color=:orange, lw=2)

    # Save graphs ----------
    savefig(mean_graph, "results/relaxed_anneal_results/$(simulation_name)_mean.png")
    savefig(mean_std_graph, "results/relaxed_anneal_results/$(simulation_name)_mean_std.png")
    savefig(specific_heat_capacity_graph, "results/relaxed_anneal_results/$(simulation_name)_specific_heat_capacity.png")
    savefig(entropy_by_temperature_graph, "results/relaxed_anneal_results/$(simulation_name)_entropy_temperature.png")
    savefig(entropy_by_absolute_energy_graph, "results/relaxed_anneal_results/$(simulation_name)_entropy_absolute.png")

    savefig(alternative_entropy_by_absolute_energy_graph, "results/relaxed_anneal_results/$(simulation_name)_alternative_entropy_absolute.png")
end


function reconstruct_histogram(histogram_data_name::String, cube)

    # Step 1: Read the CSV file into a DataFrame
    histogram_data = CSV.read(histogram_data_name, DataFrame)

    # Step 2: Extract the bin edges and weights from the DataFrame
    x_edge_starts = histogram_data.E0
    y_edge_starts = histogram_data.En
    weights_flat = histogram_data.normalised_by_E0_weights

    x_edges_unique = unique(x_edge_starts)
    y_edges_unique = unique(y_edge_starts)

    num_bins_x = length(x_edges_unique)
    num_bins_y = length(y_edges_unique)

    # Step 4: Reshape the weights into a 2D array
    weights_matrix = reshape(weights_flat, num_bins_x, num_bins_y)

    # Step 5: Construct the Histogram object with the inferred edges and reshaped weights
    # We assume the last edge for each dimension is one step beyond the last unique edge found
    # This assumes the bin width is consistent and equal to the difference between consecutive edges
    x_edges_complete = [x_edges_unique; x_edges_unique[end] + (x_edges_unique[2] - x_edges_unique[1])]
    y_edges_complete = [y_edges_unique; y_edges_unique[end] + (y_edges_unique[2] - y_edges_unique[1])]
    reconstructed_edges = (x_edges_complete, y_edges_complete)

    # Create the Histogram object
    hist = Histogram(reconstructed_edges, weights_matrix)

    return hist

end








            
