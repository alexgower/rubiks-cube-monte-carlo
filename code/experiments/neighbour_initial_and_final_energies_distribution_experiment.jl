# Xeon Only ---
# DEPOT_PATH[1]="/u/rscratch/apg59/.julia"
# using Pkg
# Pkg.instantiate()

using Plots
using GMT
using DelimitedFiles
using LaTeXStrings
include("../probes/relaxed_anneal.jl")
include("../tools/neighbour_initial_and_final_energies_graphs_plotter.jl")


# Parallel Case
# using Distributed
# addprocs(Sys.CPU_THREADS-1)
# using SharedArrays
# @everywhere include("../probes/relaxed_anneal.jl")
# @everywhere include("../tools/neighbour_initial_and_final_energies_graphs_plotter.jl")


@inbounds @fastmath function neighbour_initial_and_final_energies_distribution_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, sample_temperatures::Vector{Float64}; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0, collecting_swap_move_neighbours::Bool=false, neighbours_per_configuration_sample_size::Int64=0, average_sample_size_per_temperature::Int64=100, inherent_disorder::Bool=false, neighbour_moves_away::Int64=1)

    sample_temperatures = sort!(sample_temperatures)

    T_1 = 10.0
    T_0 = 0.09
    T_swap = 3.0
    temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
    temperature_vector = [temperature_vector; sample_temperatures]
    sort!(temperature_vector, rev=true)

    # Anneal ----------
    # Create a Rubik's cube object and run annealing function on it
    cube = RubiksCube(L)

    if inherent_disorder
        facelets = reduce(vcat, [fill(i,L^2) for i in 1:6])
        shuffle!(facelets)
        new_faces = reshape(facelets, 6, L, L)
        for i in 1:6
            cube.configuration[i][:,:] .= new_faces[i,:,:]
        end
    end

    if relaxation_iterations == 0
        temperature_vector, E_average_by_temperature, _, _, _, _, neighbour_initial_and_final_energies_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=sample_temperatures, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_initial_and_final_energies=true, average_sample_size_per_temperature=average_sample_size_per_temperature, neighbour_moves_away=neighbour_moves_away)
    else
        relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
        temperature_vector, E_average_åby_temperature, _, _, _, _, neighbour_initial_and_final_energies_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=true, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=sample_temperatures, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_initial_and_final_energies=true, average_sample_size_per_temperature=average_sample_size_per_temperature, neighbour_moves_away=neighbour_moves_away)
    end

    ##### ----- ANALYSE RESULTS -----

    printstyled("Analyzing results \n", color=:light_blue)

    # First combine (E_1, E_2) samples over all temperatures
    neighbour_initial_and_final_energies = neighbour_initial_and_final_energies_by_temperature[:]

    # Convert the array of tuples to a 2-column matrix
    data_matrix = hcat([x[1] for x in neighbour_initial_and_final_energies], [x[2] for x in neighbour_initial_and_final_energies])

    # Make energies with respect to solved energy
    data_matrix[:] .= data_matrix[:] .- solved_configuration_energy(cube)

    data_simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)"

    ### --- Save the results ---
    try
        touch(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use*"_raw_energy_connections"))

        open(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use), "w") do simulation_file
            write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Average Sample Size = $(average_sample_size_per_temperature), Neighbours Per Configuration Sample Size=$(neighbours_per_configuration_sample_size), Collecting Swap Move Neighbours = $(collecting_swap_move_neighbours) \n")
            write(simulation_file, "Initial Energy E1, Neighbour Energy E2 \n")
            
            for index in 1:size(data_matrix,1)
                write(simulation_file, "$(data_matrix[index,1]),$(data_matrix[index,2]) \n")
            end
        end

    catch ex

        println("Cannot save results to file")
        showerror(stdout, ex)

    end



    ### --- Create and Display the Graphs ---
    connectivity = collecting_swap_move_neighbours ? "Swap" : "Slice"
    neighbour_initial_and_final_energies_graph_plotter(data_simulation_name_to_use, neighbour_moves_away, connectivity)

end

# TODO DELETE BELOW IF ABOVE WORKS final_iteration_number
# ### HISTOGRAM GENERAL CALCULATIONS

#     # Determine the bin edges based on the data
#     edges = (0:1:Int(-solved_configuration_energy(cube)), 0:1:Int(-solved_configuration_energy(cube)))

#     # Compute the 2D histogram
#     hist_2d = fit(Histogram, (data_matrix[:,1], data_matrix[:,2]), edges)

#     # Normalize the histogram counts so the total volume under the histogram is 1
#     # total_volume = sum(hist_2d.weights)
#     # normalized_hist_2d = Histogram(hist_2d.edges, hist_2d.weights ./ total_volume)

#     # Normalize the histogram counts by E_1 slices
#     normalized_weights = zeros(Float64, size(hist_2d.weights))  # Create a copy to store the normalized weights
#     # Iterate through each E1 slice
#     for i in 1:size(hist_2d.weights, 1)
#         slice_sum = sum(hist_2d.weights[i, :])  # Sum of all values in the current E1 slice
#         if slice_sum != 0  # Avoid division by zero
#             normalized_weights[i, :] = hist_2d.weights[i, :] / slice_sum  # Normalize the current slice
#         end
#     end
#     normalized_hist_2d = Histogram(hist_2d.edges, normalized_weights)



#     ### --- GMT 3D HISTOGRAM ---

#     # Flatten the histogram data for GMT
#     x = repeat(collect(edges[1][1:end-1]), outer=length(edges[2])-1)
#     y = repeat(collect(edges[2][1:end-1]), inner=length(edges[1])-1)
#     z = normalized_hist_2d.weights[:]

#     ## NOTE THAT WE MAKE A LEFT HANDED SET I.E. DO A (y,x,z) PLOT TO MAKE LIKE CLAUDIO ORIGINAL PLOT
#     # Create a grid from the histogram data
#     grid_data = hcat(y, x, z)

#     # Determine the region bounds based on the data
#     region_bounds = (0, -solved_configuration_energy(cube), 0, -solved_configuration_energy(cube))

#     # Create a new grid with the normalized data
#     normalized_grid = xyz2grd(grid_data, region=region_bounds, inc=1)

#     # Plot the normalized 3D histogram and display it immediately
#     bar3(normalized_grid, view=(190,50), frame=(xlabel="E2", ylabel="E1", zlabel="Frequency", annot=:auto, ticks=:auto, grid=:auto), fill=[0,115,190], lw=0.25, fmt=:png, show=true, savefig="results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_Ps=$(swap_move_probability)_3D.png")



#     ### --- Plots.jl 2D HISTOGRAM ---

#     # Create the plot using histogram2d
#     bin_edges_x = 0:1:-solved_configuration_energy(cube)
#     bin_edges_y = 0:1:-solved_configuration_energy(cube)

#     # Compute the 2D histogram
#     hist_2d = fit(Histogram, (data_matrix[:,1], data_matrix[:,2]), edges)

#     # Initialize an array to hold the biases
#     biases = zeros(size(data_matrix, 1))

#     # NOTE THAT BIASES BELOW DO NOT MEAN FREQUENCIES IN THE HISTOGRAM BUT RATHER HOW WE BIAS THESE FREQUENCIES
#     # Calculate the sum of counts in each E_1 slice
#     for i in 1:size(hist_2d.weights, 1)
#         slice_sum = sum(hist_2d.weights[i, :])
#         if slice_sum != 0  # Avoid division by zero
#             # Assign the reciprocal of the slice sum to the weights array for each data point in the current E_1 slice
#             mask = data_matrix[:,1] .== bin_edges_x[i]
#             biases[mask] .= 1.0 / slice_sum
#         end
#     end

#     ## NOTE THAT WE AGAIN PLOT E_2 ON THE X-AXIS FOR CONSISTENCY

#     graph = histogram2d(data_matrix[:,2], data_matrix[:,1], show_empty_bins=true,
#     normalize=:pdf, bins=(bin_edges_y, bin_edges_x), weights=biases, color=:plasma, xlabel="E2", ylabel="E1", xlims=(0, -solved_configuration_energy(cube)), ylims=(0, -solved_configuration_energy(cube)), zlabel="Frequency")



#     ### --- Scatter Plot of Modal E_2 for each E_1 ---

#     ## GET MODAL POINTS

#     # Get the unique E1 values
#     E1_bin_values = bin_edges_x[1:end-1]  # Exclude the last bin edge as it does not correspond to a bin center
#     modal_E2 = zeros(Float64, length(E1_bin_values))

#     # Iterate through each E1 slice to find the modal E2 value
#     for (i, E1) in pairs(E1_bin_values)
#         # Get the weights for the current E1 slice
#         weights_slice = hist_2d.weights[i, :]
        
#         # Find the index of the modal E2 value in the current E1 slice
#         modal_index = argmax(weights_slice)
#         # Convert the index to the modal E2 value using the bin edges
#         modal_value = bin_edges_y[modal_index]
#         # Store the modal E2 value
#         modal_E2[i] = modal_value

#     end 

#     # Filter out all modal_E2 values that are 0.0 
#     filtered_indices = findall(modal_E2 .!= 0.0)
#     modal_E1_values = E1_bin_values[filtered_indices]
#     modal_E2 = modal_E2[filtered_indices]




#     ## GET LOWER TAIL POINTS

#     E1_values = unique(data_matrix[:,1])
#     lower_tail_E2 = zeros(Float64, length(E1_values))

#     # Iterate through each E1 slice to find the modal E2 value
#     for (i, E1) in pairs(E1_values)

#         # Get the weights for the current E1 slice
#         E_2_values = [data_matrix[j,2] for j in 1:size(data_matrix,1) if data_matrix[j,1] == E1]

#         # Find the index of the lower tail E2 value in the current E1 slice
#         lower_tail = minimum(E_2_values)
#         lower_tail_E2[i] = lower_tail

#     end 

#     ##  CREATE THE SCATTER PLOT
#     mode_graph = Plots.scatter(lower_tail_E2, E1_values, xlabel="E2", ylabel="E1", label="Minimum E2", color=:red)
#     Plots.scatter!(mode_graph, modal_E2, modal_E1_values, label="Modal E2", color=:blue)
#     # legend!(mode_graph, :bottomright)


#     # Add E_1 = E_2 lines to graphs
#     min_value = minimum([minimum(data_matrix[:,1]), minimum(data_matrix[:,2])])
#     max_value = maximum([maximum(data_matrix[:,1]), maximum(data_matrix[:,2])])
#     Plots.plot!(graph, [0, max_value], [0, max_value], line=:dash, color=:green, lw=2, legend=false)
#     Plots.plot!(mode_graph, [0, max_value], [0, max_value], line=:dash, color=:green, lw=2, legend=false)

#     # Add E_c lines to graphs (average energy at onset of transition i.e. E(T_c^+))
#     E_c = 95.0
#     if min_value < E_c < max_value
#         vline!(graph, [E_c], line=:dash, color=:black, lw=1)
#         vline!(mode_graph, [E_c], line=:dash, color=:black, lw=1)
#     end


#     ### --- Save and display the graphs ---
#     savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_Ps=$(swap_move_probability)_2D.png")
#     savefig(mode_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(simulation_name)_Ps=$(swap_move_probability)_2D_mode.png")



















######## NOT FROM HEREE


# @everywhere @inbounds @fastmath function parallel_neighbour_energy_deltas_distribution_experiment(simulation_name::String, L::Int64, swap_move_probability::Float64, N_T::Int64, neighbour_sample_temperatures::Vector{Float64}, average_sample_size_per_temperature::Int64; verbose_metropolis_swap::Bool=false, normalization::String="solved", relaxation_iterations::Int64=Int(0), mixing_p_swap::Float64=0.0, collecting_swap_move_neighbours::Bool=false, neighbours_per_configuration_sample_size::Int64=0, collect_minimum_neighbour_energy_info_only::Bool=false, extra_swap_moves::Int64=0, extra_slice_rotations::Int64=0)

#     neighbour_sample_temperatures = sort!(neighbour_sample_temperatures)

#     T_1 = 10.0
#     T_0 = 0.09
#     T_swap = 1.5
#     temperature_vector::Vector{Float64} = [T_1*(T_0/T_1)^(m/N_T) for m in 0:N_T]
#     temperature_vector = [temperature_vector; neighbour_sample_temperatures]
#     sort!(temperature_vector, rev=true)

#     if neighbours_per_configuration_sample_size == 0
#         neighbours_per_configuration_sample_size = configuration_network_degree(L, false)
#     end

#     # Anneal ----------
#     # Create a Rubik's cube object and run annealing function on it
#     cube = RubiksCube(L)

#     E_values_by_temperature = SharedArray{Float64,2}(length(temperature_vector), average_sample_size_per_temperature)

#     if !collect_minimum_neighbour_energy_info_only
#         neighbour_energy_deltas_by_temperature = SharedArray{Float64,2}(length(neighbour_sample_temperatures), neighbours_per_configuration_sample_size*average_sample_size_per_temperature)
#     else
#         neighbour_energy_deltas_by_temperature = SharedArray{Float64,2}(length(neighbour_sample_temperatures), 1*average_sample_size_per_temperature)
#     end

#     # neighbour_energy_deltas_by_temperature .= ones(length(neighbour_sample_temperatures), neighbours_per_configuration_sample_size*average_sample_size_per_temperature)
#     # neighbour_energy_deltas_by_temperature = zeros(length(neighbour_sample_temperatures), neighbours_per_configuration_sample_size*average_sample_size_per_temperature)


#     @sync @distributed for trial in 1:average_sample_size_per_temperature
#         printstyled("Trial: $trial \n", color=:light_blue)
#         if relaxation_iterations == 0
#             temperature_vector, trial_E_average_by_temperature, _, _, _, _, trial_neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=false, verbose_metropolis_swap=verbose_metropolis_swap, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=neighbour_sample_temperatures, average_sample_size_per_temperature=1, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_minimum_neighbour_energy_info_only=collect_minimum_neighbour_energy_info_only, extra_swap_moves=extra_swap_moves, extra_slice_rotations=extra_slice_rotations)
#         else
#             relaxation_iterations_vector = [relaxation_iterations for T in temperature_vector]
#             temperature_vector, trial_E_average_by_temperature, _, _, _, _, trial_neighbour_energy_deltas_by_temperature = relaxed_anneal!(cube, temperature_vector; swap_move_probability=swap_move_probability, T_swap=T_swap, verbose_annealing=false, verbose_metropolis_swap=verbose_metropolis_swap, relaxation_iterations_vector = relaxation_iterations_vector, mixing_p_swap=mixing_p_swap, neighbour_sample_temperatures=neighbour_sample_temperatures, average_sample_size_per_temperature=1, collecting_swap_move_neighbours=collecting_swap_move_neighbours, neighbours_per_configuration_sample_size=neighbours_per_configuration_sample_size, collect_minimum_neighbour_energy_info_only=collect_minimum_neighbour_energy_info_only, extra_swap_moves=extra_swap_moves, extra_slice_rotations=extra_slice_rotations)
#         end

#         if !collect_minimum_neighbour_energy_info_only
#             neighbour_energy_deltas_by_temperature[:,(trial-1)*neighbours_per_configuration_sample_size+1:trial*neighbours_per_configuration_sample_size] .= trial_neighbour_energy_deltas_by_temperature[:,:]
#         else
#             neighbour_energy_deltas_by_temperature[:,trial] .= trial_neighbour_energy_deltas_by_temperature[:,:]
#         end
#         E_values_by_temperature[:,trial] .= trial_E_average_by_temperature[:]
#     end

#     E_average_by_temperature = [mean(E_values_by_temperature[index,:]) for index in 1:length(temperature_vector)]


#     # Analyse Results ---------

#     # Find maximum and minima neighbour energy deltas to use as end points of histogram
#     max_neighbour_energy_delta = maximum(neighbour_energy_deltas_by_temperature)
#     min_neighbour_energy_delta = minimum(neighbour_energy_deltas_by_temperature)
    
#     # Find maximimum probability density
#     max_probability_density = 0
#     for (index,sample_temperature) in pairs(neighbour_sample_temperatures)
#         histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])
#         probability_densities = collect(values(histogram))
#         probability_densities = probability_densities ./ sum(probability_densities)
#         this_max_probability_density = max(max_probability_density, maximum(probability_densities))

#         if this_max_probability_density > max_probability_density
#             max_probability_density = this_max_probability_density
#         end
#     end

#     for (index,sample_temperature) in pairs(neighbour_sample_temperatures)

#         simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)_T=$(sample_temperature)"
        
#         # Create plot ----------

#         # Find maximum and minima neighbour energy deltas to use as end points of histogram
#         max_neighbour_energy_delta = maximum(neighbour_energy_deltas_by_temperature)
#         min_neighbour_energy_delta = minimum(neighbour_energy_deltas_by_temperature)

#         # Find maximimum probability density
#         max_probability_density = 0
#         for (index,sample_temperature) in pairs(neighbour_sample_temperatures)
#             histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])
#             probability_densities = collect(values(histogram))
#             probability_densities = probability_densities ./ sum(probability_densities)
#             this_max_probability_density = max(max_probability_density, maximum(probability_densities))

#             if this_max_probability_density > max_probability_density
#                 max_probability_density = this_max_probability_density
#             end
#         end

#         try
            
#             # Create a histogram
#             histogram = countmap(neighbour_energy_deltas_by_temperature[index,:])

#             # Now add a bunch of zeros to the histogram for all range of energy deltas between min_neighbour_energy_delta and max_neighbour_energy_delta so that all histograms have same width
#             for energy_delta in min_neighbour_energy_delta:max_neighbour_energy_delta
#                 if !haskey(histogram, energy_delta)
#                     histogram[energy_delta] = 0
#                 end
#             end

#             energy_deltas = collect(keys(histogram))
#             println("Energy deltas: ", energy_deltas)

#             probability_densities = collect(values(histogram))
#             println("Probability densities: ", probability_densities)

#             # Normalize the values
#             probability_densities = probability_densities ./ sum(probability_densities)

#             # Plot histogram
#             title = collect_minimum_neighbour_energy_info_only ? "Min Neighbour Energy Deltas Distribution at T=$(sample_temperature)" : "Neighbour Energy Deltas Distribution at T=$(sample_temperature)"
#             label = collect_minimum_neighbour_energy_info_only ? "Min Neighbour Energy Deltas" : "Neighbouring Energy Deltas"
#             graph = bar(energy_deltas, probability_densities, xticks=energy_deltas, xtickfontsize=4, label=label, ylim=(0,max_probability_density*1.1), legend=:topleft, title=title, xlabel=label, ylabel="Probability Density")

#             # Plot One in a Hundred Energy Dashed Line
#             one_in_hundred_energy = log(100)*neighbour_sample_temperatures[index]
#             vline!(graph, [one_in_hundred_energy], linestyle=:dash, label="1/100 Energy", color=:orange)

#             # Plot One in a Thousand Energy Dashed Line
#             one_in_thousand_energy = log(1000)*neighbour_sample_temperatures[index]
#             vline!(graph, [one_in_thousand_energy], linestyle=:dash, label="1/1000 Energy", color=:red)

#             # Save graphs ----------
#             savefig(graph, "results/neighbour_energy_deltas_distribution_results/$(simulation_name_to_use).png")


#         catch ex

#             println("Cannot display or save plot")
#             showerror(stdout, ex)

#         end
#     end

#     # Save Neighbour Energy Delta Results ----------
#     try
#         data_simulation_name_to_use = "$(simulation_name)_Ps=$(swap_move_probability)"

#         touch(joinpath("results/neighbour_energy_deltas_distribution_results",data_simulation_name_to_use))

#         open(joinpath("results/neighbour_energy_deltas_distribution_results",data_simulation_name_to_use), "w") do simulation_file
#             write(simulation_file, "Simulation:L=$L, P_s=$swap_move_probability, T_swap=$T_swap, T_1=$T_1, T_0=$T_0, N_T=$N_T, Trials=$(average_sample_size_per_temperature), Neighbours Per Configuration Sample Size=$(neighbours_per_configuration_sample_size), Collecting Swap Move Neighbours = $(collecting_swap_move_neighbours), Average Sample Size = $(average_sample_size_per_temperature) \n")
#             write(simulation_file, "Sample Temperature T, Average Energy E, Neighbour Energy Deltas \n")
            
#             for (index,sample_temperature) in pairs(neighbour_sample_temperatures)
#                 write(simulation_file, "$(neighbour_sample_temperatures[index]),$(E_average_by_temperature[index]), $(neighbour_energy_deltas_by_temperature[index,:]) \n")
#             end
#         end

#     catch ex

#         println("Cannot save results to file")
#         showerror(stdout, ex)

#     end



# end






