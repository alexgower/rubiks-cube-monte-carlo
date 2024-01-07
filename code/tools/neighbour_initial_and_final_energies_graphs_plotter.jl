using DelimitedFiles
using Plots

using StatsBase
using LaTeXStrings
using CSV
using DataFrames

using Colors
using ColorTypes, ColorSchemes
using Plots.PlotMeasures

include("../core/rubiks_cube.jl")

function neighbour_initial_and_final_energies_graphs_plotter(data_simulation_name_to_use::String, connectivity::String="Slice"; neighbour_moves_away::Int64=1, E_star::Float64=0.0)

    n = neighbour_moves_away

    ### -- READ IN THE DATA --
    header_line = readlines(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use))[1]
    data_matrix = readdlm(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use), ',', Float64, '\n', skipstart=2)

    ### -- SET UP DEFAULT PARAMETERS --
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    # Remove NaNs and overflow # TODO how exist?
    data_matrix = remove_bad_rows(data_matrix, L)

    ### -- HISTOGRAM GENERAL CALCULATIONS --

    # Determine the bin edges based on the data
    edges = (0:1:Int(-solved_configuration_energy(cube)), 0:1:Int(-solved_configuration_energy(cube)))

    # Compute the 2D histogram
    hist_2d = fit(Histogram, (data_matrix[:,1], data_matrix[:,2]), edges)

    # Normalize the histogram counts so the total volume under the histogram is 1
    # total_volume = sum(hist_2d.weights)
    # normalized_hist_2d = Histogram(hist_2d.edges, hist_2d.weights ./ total_volume)

    # Normalize the histogram counts by E_0 slices
    normalized_weights = zeros(Float64, size(hist_2d.weights))  # Create a copy to store the normalized weights
    # Iterate through each E0 slice
    for i in 1:size(hist_2d.weights, 1)
        slice_sum = sum(hist_2d.weights[i, :])  # Sum of all values in the current E0 slice
        println("Number of weights at a single E_0 = $i ", slice_sum)
        if slice_sum != 0  # Avoid division by zero
            normalized_weights[i, :] = hist_2d.weights[i, :] / slice_sum  # Normalize the current slice
        end
    end

    normalized_hist_2d = Histogram(hist_2d.edges, normalized_weights)








    ### -- SAVE HISTOGRAM DATA --

    # Flatten the histogram weights
    # Assume hist_2d is your histogram object
    # Flatten the 2D weights array into a 1D array
    weights_flat = vec(hist_2d.weights)
    normalised_weights_flat = vec(normalized_hist_2d.weights)

    # Generate all combinations of x and y bin edges that correspond to the weights
    # Note that we exclude the last edge since it is not the start of a bin
    x_bin_edges = repeat(hist_2d.edges[1][1:end-1], outer=length(hist_2d.edges[2])-1)
    y_bin_edges = repeat(hist_2d.edges[2][1:end-1], inner=length(hist_2d.edges[1])-1)

    # Create a DataFrame with the correct edge pairing and flattened weights
    histogram_data = DataFrame(
        E0 = x_bin_edges,
        E_n = y_bin_edges,
        weights = weights_flat
    )


    # Create a DataFrame with the correct edge pairing and flattened weights
    normalised_histogram_data = DataFrame(
        E0 = x_bin_edges,
        En = y_bin_edges,
        normalised_by_E0_weights = normalised_weights_flat
    )

    # Save the DataFrames to CSV files
    CSV.write(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use*"_histogram_data.csv"), histogram_data)
    CSV.write(joinpath("results/neighbour_initial_and_final_energies_distribution_results",data_simulation_name_to_use*"_normalised_histogram_data.csv"), normalised_histogram_data)






    # ### -- GMT 3D HISTOGRAM --

    # # Flatten the histogram data for GMT
    # x = repeat(collect(edges[1][1:end-1]), outer=length(edges[2])-1)
    # y = repeat(collect(edges[2][1:end-1]), inner=length(edges[1])-1)
    # z = normalized_hist_2d.weights[:]

    # ## NOTE THAT WE MAKE A LEFT HANDED SET I.E. DO A (y,x,z) PLOT TO MAKE LIKE CLAUDIO ORIGINAL PLOT
    # # Create a grid from the histogram data
    # grid_data = hcat(y, x, z)

    # # Determine the region bounds based on the data
    # region_bounds = (0, -solved_configuration_energy(cube), 0, -solved_configuration_energy(cube))

    # # Create a new grid with the normalized data
    # normalized_grid = xyz2grd(grid_data, region=region_bounds, inc=1)

    # # Plot the normalized 3D histogram and display it immediately
    # bar3(normalized_grid, view=(190,50), FONT="12p,Times-Roman", FONT_TITLE="10p", MAP_TITLE_OFFSET="0/2c", frame=(xlabel="E@-n@-", ylabel="E@-0@-", zlabel="Frequency", annot=:auto, ticks=:auto, grid=:auto, title= " (E@-0@-, "*"E@-n@-)"* " (n=$n) $connectivity Connectivity Frequency Histogram L=$L"), fill=[0,115,190], lw=0.25, fmt=:png, show=true,  savefig="results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_3D.png")










    min_value = minimum([minimum(data_matrix[:,1]), minimum(data_matrix[:,2])])
    max_value = maximum([maximum(data_matrix[:,1]), maximum(data_matrix[:,2])])

    ### -- Plots. jl 2D HISTOGRAM --

    # Create the plot using histogram2d
    bin_edges_x = 0:1:-solved_configuration_energy(cube)
    bin_edges_y = 0:1:-solved_configuration_energy(cube)

    # Compute the 2D histogram
    hist_2d = fit(Histogram, (data_matrix[:,1], data_matrix[:,2]), edges)

    # Initialize an array to hold the biases
    biases = zeros(size(data_matrix, 1))

    # NOTE THAT BIASES BELOW DO NOT MEAN FREQUENCIES IN THE HISTOGRAM BUT RATHER HOW WE BIAS THESE FREQUENCIES
    # Calculate the sum of counts in each E_0 slice
    for i in 1:size(hist_2d.weights, 1)
        slice_sum = sum(hist_2d.weights[i, :])
        if slice_sum != 0  # Avoid division by zero
            # Assign the reciprocal of the slice sum to the weights array for each data point in the current E_0 slice
            mask = data_matrix[:,1] .== bin_edges_x[i]
            biases[mask] .= 1.0 / slice_sum
        end
    end


    graph = histogram2d(
    data_matrix[:,1]./-solved_configuration_energy(cube), 
    data_matrix[:,2]./-solved_configuration_energy(cube), 
    color=:bluesreds, 
    show_empty_bins=false,
    normalize=:pdf, 
    bins=(bin_edges_x./-solved_configuration_energy(cube), bin_edges_y./-solved_configuration_energy(cube)), 
    weights=biases, 
    ylabel="Neighbour Configuration Energy, "*L"E^{(1)}", 
    xlabel="Initial Configuration Energy, "*L"E^{(0)}", 
    title="L=$L $connectivity Cube Energy Connectivity", 
    xlims=(min_value./-solved_configuration_energy(cube), max_value./-solved_configuration_energy(cube)), 
    ylims=(min_value./-solved_configuration_energy(cube), max_value./-solved_configuration_energy(cube)), 
    colorbar_title="Sampled Frequency",
    titlefontsize=10,   # Title font size
    xguidefontsize=8,   # X-axis label font size
    yguidefontsize=8,   # Y-axis label font size
    margin=5mm          # Margin around the plot
)
    
    


    ### -- Scatter Plot of Modal E_n for each E_0 --

    ## GET MODAL POINTS

    # Get the unique E0 values
    E0_bin_values = bin_edges_x[1:end-1]  # Exclude the last bin edge as it does not correspond to a bin center
    modal_En = zeros(Float64, length(E0_bin_values))

    # Iterate through each E0 slice to find the modal En value
    for (i, E0) in pairs(E0_bin_values)
        # Get the weights for the current E0 slice
        weights_slice = hist_2d.weights[i, :]
        
        # Find the index of the modal En value in the current E0 slice
        modal_index = argmax(weights_slice)
        # Convert the index to the modal En value using the bin edges
        modal_value = bin_edges_y[modal_index]
        # Store the modal En value
        modal_En[i] = modal_value

    end 

    # Filter out all modal_En values that are 0.0 
    filtered_indices = findall(modal_En .!= 0.0)
    modal_E0_values = E0_bin_values[filtered_indices]
    modal_En = modal_En[filtered_indices]




    ## GET LOWER AND UPPER TAIL POINTS

    E0_values = unique(data_matrix[:,1])
    lower_tail_En = zeros(Float64, length(E0_values))
    upper_tail_En = zeros(Float64, length(E0_values))

    # Iterate through each E0 slice to find the modal En value
    for (i, E0) in pairs(E0_values)

        # Get the weights for the current E0 slice
        E_n_values = [data_matrix[j,2] for j in 1:size(data_matrix,1) if data_matrix[j,1] == E0] 

        # Find the index of the lower tail En value in the current E0 slice
        lower_tail = minimum(E_n_values)
        lower_tail_En[i] = lower_tail

        # Find the index of the upper tail En value in the current E0 slice
        upper_tail = maximum(E_n_values)
        upper_tail_En[i] = upper_tail

    end 

 
    lower_tail_En_below_E0 = lower_tail_En[lower_tail_En .< E0_values]
    E_0_values_for_lower_tail_En_below_E0 = E0_values[lower_tail_En .< E0_values]

    lower_tail_En_equal_E0 = lower_tail_En[lower_tail_En .== E0_values]
    E_0_values_for_lower_tail_En_equal_E0 = E0_values[lower_tail_En .== E0_values]

    lower_tail_En_above_E0 = lower_tail_En[lower_tail_En .> E0_values]
    E_0_values_for_lower_tail_En_above_E0 = E0_values[lower_tail_En .> E0_values]

    ##  CREATE THE SCATTER PLOT
    mode_graph = plot(
        ylabel="Neighbour Configuration Energy, "*L"E^{(1)}", 
        xlabel="Initial Configuration Energy, "*L"E^{(0)}", 
        title="L=$L $connectivity Cube Energy Connectivity", 
        legend=:bottomright, 
        legendfont = font(8, "Times"), 
        titlefontsize=10,   # Title font size
        xguidefontsize=8,   # X-axis label font size
        yguidefontsize=8,   # Y-axis label font size
        margin=5mm          # Margin around the plot
    )    
    scatter!(mode_graph, E_0_values_for_lower_tail_En_below_E0./-solved_configuration_energy(cube), lower_tail_En_below_E0./-solved_configuration_energy(cube), label="Minimal E⁽¹⁾<E⁽⁰⁾ for E⁽⁰⁾", color=:red)
    scatter!(mode_graph, E_0_values_for_lower_tail_En_equal_E0./-solved_configuration_energy(cube), lower_tail_En_equal_E0./-solved_configuration_energy(cube),  label="Minimal E⁽¹⁾=E⁽⁰⁾ for E⁽⁰⁾", color=RGB(0.6, 1.0, 0.6))
    scatter!(mode_graph, E_0_values_for_lower_tail_En_above_E0./-solved_configuration_energy(cube), lower_tail_En_above_E0./-solved_configuration_energy(cube), label="Minimal E⁽¹⁾>E⁽⁰⁾ for E⁽⁰⁾", color=:orange)
    scatter!(mode_graph, modal_E0_values./-solved_configuration_energy(cube), modal_En./-solved_configuration_energy(cube), label="Modal E⁽¹⁾ for E⁽⁰⁾", color=:blue)

    lower_tail_graph = plot(
        ylabel="Neighbour Configuration Energy, "*L"E^{(1)}", 
        xlabel="Initial Configuration Energy, "*L"E^{(0)}", 
        title="L=$L $connectivity Cube Energy Connectivity", 
        legend=:bottomright, 
        legendfont = font(8, "Times"), 
        titlefontsize=10,   # Title font size
        xguidefontsize=8,   # X-axis label font size
        yguidefontsize=8,   # Y-axis label font size
        margin=5mm          # Margin around the plot
    )
    scatter!(lower_tail_graph, E_0_values_for_lower_tail_En_below_E0./-solved_configuration_energy(cube), lower_tail_En_below_E0./-solved_configuration_energy(cube), label="Minimal E⁽¹⁾<E⁽⁰⁾ for E⁽⁰⁾", color=:red)
    scatter!(lower_tail_graph, E_0_values_for_lower_tail_En_equal_E0./-solved_configuration_energy(cube), lower_tail_En_equal_E0./-solved_configuration_energy(cube),  label="Minimal E⁽¹⁾=E⁽⁰⁾ for E⁽⁰⁾", color=RGB(0.6, 1.0, 0.6))
    scatter!(lower_tail_graph, E_0_values_for_lower_tail_En_above_E0./-solved_configuration_energy(cube), lower_tail_En_above_E0./-solved_configuration_energy(cube), label="Minimal E⁽¹⁾>E⁽⁰⁾ for E⁽⁰⁾", color=:orange)

    # Add E_0 = E_n lines to graphs
    plot!(graph, [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], line=:dash, color=:orange, lw=2, label="E⁽⁰⁾=E⁽¹⁾")
    plot!(mode_graph, [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], line=:dash, color=:orange, lw=2, label="E⁽⁰⁾=E⁽¹⁾")
    plot!(lower_tail_graph, [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], line=:dash, color=:black, lw=2, label="E⁽⁰⁾=E⁽¹⁾")

    # Add E_star vertical lines to graphs if E_star is nonzero
    if E_star != 0.0
        vline!(graph, [E_star], line=:dash, color=:green, lw=2, label="")
        vline!(mode_graph, [E_star], line=:dash, color=:green, lw=2, label="")
        vline!(lower_tail_graph, [E_star], line=:dash, color=:green, lw=2, label="")
        annotate!(graph, [(E_star+0.02, ylims(graph)[1]+0.02, Plots.text(L"E^*", 8, :black))])
        annotate!(mode_graph, [(E_star+0.02, ylims(mode_graph)[1]+0.02, Plots.text(L"E^*", 8, :black))])
        annotate!(lower_tail_graph, [(E_star+0.02, ylims(lower_tail_graph)[1]+0.02, Plots.text(L"E^*", 8, :black))])
    end

    ### -- Save and display the graphs --
    println("Saving graphs...")
    savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D.svg")
    savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D.png")
    display(graph)
    println("Saved histogram")

    savefig(mode_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D_mode.svg")
    savefig(mode_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D_mode.png")
    display(mode_graph)

    # # and maximal graph tools
    scatter!(mode_graph, E0_values./-solved_configuration_energy(cube), upper_tail_En./-solved_configuration_energy(cube),  label="Maximum E⁽¹⁾ for E⁽⁰⁾", color=:green)
    savefig(mode_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D_mode_maximal.svg")
    savefig(mode_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D_mode_maximal.png")
    display(mode_graph)

    # and lower tail graph
    savefig(lower_tail_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D_lower_tail.svg")
    savefig(lower_tail_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(data_simulation_name_to_use)_2D_lower_tail.png")
    display(lower_tail_graph)

    println("Saved everything")
    println("E_0 values where Lower Tail E_1 is ABOVE E_0: ", sort(E_0_values_for_lower_tail_En_above_E0))
    println("E_0 values where Lower Tail E_1 is EQUAL TO E_0: ", sort(E_0_values_for_lower_tail_En_equal_E0))

end

function remove_bad_rows(data::Array{Float64,2}, L::Int64)::Array{Float64,2}
    # Find rows without NaN or negative values or above solved configuration energy in the first two columns or above solved configuration energy
    non_bad_rows = .!isnan.(data[:, 1]) .& .!isnan.(data[:, 2]) .& (data[:, 1] .>= 0) .& (data[:, 2] .>= 0) .& (data[:, 1] .< -solved_configuration_energy(RubiksCube(L))) .& (data[:, 2] .< -solved_configuration_energy(RubiksCube(L)))
    # Return the data without rows containing NaN
    return data[non_bad_rows, :]
end