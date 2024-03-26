using DelimitedFiles
using LaTeXStrings, Plots; 

using StatsBase
using CSV
# using DataFrames

using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

include("../core/rubiks_cube.jl")

# Make nice function renamings
E1_E2_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation") = energy_connectivity_histogram_figure(simulation_name, connectivity; neighbour_order_to_measure_to=2)
E0_E1_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation") = energy_connectivity_histogram_figure(simulation_name, connectivity; neighbour_order_to_measure_to=1)



# Main function
function energy_connectivity_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation"; neighbour_order_to_measure_to::Int64=1)

    E_star = -0.39015151515151514

    ### --- READ IN THE DATA ---
    filename = "results/final_paper_results",simulation_name*"_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_connections"

    header_line = readlines(joinpath(filename))[1]
    energy_connections_data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=2)


    ### --- SET UP DEFAULT PARAMETERS ---

    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    # Remove NaNs and overflow # TODO how exist?
    energy_connections_data_matrix = remove_bad_rows(energy_connections_data_matrix, L)
    Plots.default(dpi = 300)





    ### --- HISTOGRAM GENERAL CALCULATIONS ---

    ## -- DIAGONAL GRAPH --

    # Determine the bin edges based on the data
    bin_edges_x = 0:1:-solved_configuration_energy(cube)
    bin_edges_y = 0:1:-solved_configuration_energy(cube)
    edges = (bin_edges_x, bin_edges_y)

    # Compute the 2D histogram
    hist_2d = fit(Histogram, (energy_connections_data_matrix[:,1], energy_connections_data_matrix[:,2]), edges)


    ### -- Plots. jl 2D HISTOGRAM --

    min_value = minimum([minimum(energy_connections_data_matrix[:,1]), minimum(energy_connections_data_matrix[:,2])])
    max_value = maximum([maximum(energy_connections_data_matrix[:,1]), maximum(energy_connections_data_matrix[:,2])])


    # Initialize an array to hold the biases
    biases = zeros(size(energy_connections_data_matrix, 1))

    # NOTE THAT BIASES BELOW DO NOT MEAN FREQUENCIES IN THE HISTOGRAM BUT RATHER HOW WE BIAS THESE FREQUENCIES
    # Calculate the sum of counts in each E_0 slice
    for i in 1:size(hist_2d.weights, 1)
        slice_sum = sum(hist_2d.weights[i, :])
        if slice_sum != 0  # Avoid division by zero
            # Assign the reciprocal of the slice sum to the weights array for each data point in the current E_0 slice
            mask = energy_connections_data_matrix[:,1] .== bin_edges_x[i]
            biases[mask] .= 1.0 / slice_sum
        end
    end

    E0_values = energy_connections_data_matrix[:,1]./-solved_configuration_energy(cube)
    E1_values = energy_connections_data_matrix[:,2]./-solved_configuration_energy(cube)

    if neighbour_order_to_measure_to == 2
        diagonal_ylabel = "Second-Neighbour Configuration Energy Density, "*L"E^{(2)}/E_s"
        diagonal_xlabel = "Neighbour Configuration Energy Density, "*L"E^{(1)}/E_s"
    elseif neighbour_order_to_measure_to == 1
        diagonal_ylabel = "Neighbour Configuration Energy Density, "*L"E^{(1)}/E_s"
        diagonal_xlabel = "Initial Configuration Energy Density, "*L"E^{(0)}/E_s"
    else
        error("Invalid neighbour_order_to_measure_to value for x/y labels: $neighbour_order_to_measure_to")
    end

    graph = histogram2d(
        E0_values, 
        E1_values, 
        color=:bluesreds, 
        show_empty_bins=false,
        normalize=:pdf, 
        bins=(bin_edges_x./-solved_configuration_energy(cube), bin_edges_y./-solved_configuration_energy(cube)), 
        weights=biases, 
        ylabel=diagonal_ylabel, 
        xlabel=diagonal_xlabel, 
        title="$connectivity Cube", 
        xlims=(minimum(E0_values), maximum(E0_values)), 
        ylims=(minimum(E1_values), maximum(E1_values)),
        colorbar_title="Sampled Frequency",
        titlefontsize=10,   # Title font size
        xguidefontsize=8,   # X-axis label font size
        yguidefontsize=8,   # Y-axis label font size
        margin=5mm          # Margin around the plot
    )

    # Add E_0 = E1 lines to graph
    if neighbour_order_to_measure_to==2
        annotation = "E⁽¹⁾=E⁽²⁾"
    elseif neighbour_order_to_measure_to==1
        annotation = "E⁽⁰⁾=E⁽¹⁾"
    else
        error("Invalid neighbour_order_to_measure_to value for annotation: $neighbour_order_to_measure_to")
    end

    plot!(graph, [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], line=:dash, color=:orange, lw=2, label=annotation)

    # Add E_star vertical lines to graphs if E_star is nonzero and we are slice rotation cube
    if E_star != 0.0 && connectivity == "Slice-Rotation"
        vline!(graph, [E_star], line=:dash, color=:green, lw=2, label="")
        annotate!(graph, [(E_star+0.02, ylims(graph)[1]+0.02, Plots.text(L"E^*", 8, :black))])
    end

    ### -- Save and display the graphs --
    println("Saving diagonal graph...")
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal.svg")
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal.png")
    display(graph)



    ## -- HORIZONTAL GRAPH --

    # Determine the bin edges based on the data
    bin_edges_x = 0:1:-solved_configuration_energy(cube)
    bin_edges_y = solved_configuration_energy(cube):1:-solved_configuration_energy(cube)
    edges = (bin_edges_x, bin_edges_y)

    # Compute the 2D histogram
    hist_2d = fit(Histogram, (energy_connections_data_matrix[:,1], energy_connections_data_matrix[:,2].-energy_connections_data_matrix[:,1]), edges)


    ### -- Plots. jl 2D HISTOGRAM --

    # Initialize an array to hold the biases
    biases = zeros(size(energy_connections_data_matrix, 1))

    # NOTE THAT BIASES BELOW DO NOT MEAN FREQUENCIES IN THE HISTOGRAM BUT RATHER HOW WE BIAS THESE FREQUENCIES
    # Calculate the sum of counts in each E_0 slice
    for i in 1:size(hist_2d.weights, 1)
        slice_sum = sum(hist_2d.weights[i, :])
        if slice_sum != 0  # Avoid division by zero
            # Assign the reciprocal of the slice sum to the weights array for each data point in the current E_0 slice
            mask = energy_connections_data_matrix[:,1] .== bin_edges_x[i]
            biases[mask] .= 1.0 / slice_sum
        end
    end

    E0_values = (energy_connections_data_matrix[:,1]./-solved_configuration_energy(cube))
    E_difference_values = (energy_connections_data_matrix[:,2] .- energy_connections_data_matrix[:,1])./-solved_configuration_energy(cube)

    if neighbour_order_to_measure_to == 2
        horizontal_ylabel = "Second-Neighbour Energy Density Difference, "*L"(E^{(2)}\!\!\!-\!E^{(1)})/\,|\!\!E_s\,|"
        horizontal_xlabel = "Neighbour Energy Density, "*L"E^{(1)}\!\!/\,|\!\!E_s\,|"
    elseif neighbour_order_to_measure_to == 1
        horizontal_ylabel = "Neighbour Energy Density Difference, "*L"(E^{(1)}\!\!\!-\!E^{(0)})/\,|\!\!E_s\,|"
        horizontal_xlabel = "Energy Density, "*L"E^{(0)}\!\!/\,|\!\!E_s\,|"
    else
        error("Invalid neighbour_order_to_measure_to value for x/y labels: $neighbour_order_to_measure_to")
    end

    graph = histogram2d(
        E0_values, 
        E_difference_values, 
        color=:bluesreds, 
        show_empty_bins=false,
        normalize=:pdf, 
        bins=((bin_edges_x./-solved_configuration_energy(cube)), bin_edges_y./-solved_configuration_energy(cube)), 
        weights=biases, 
        ylabel=horizontal_ylabel, 
        xlabel=horizontal_xlabel,
        xlims=(minimum(E0_values), maximum(E0_values)), 
        ylims=(minimum(E_difference_values), maximum(E_difference_values)),
        colorbar_title="Sampled Frequency",
        xguidefontsize=8,   # X-axis label font size
        yguidefontsize=8,   # Y-axis label font size
        margin=1mm,          # Margin around the plot
        title="$connectivity Cube",
        titlefontsize=10,   # Title font size
        xticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
    )


    # Add E_0 = E1 lines to graph
    hline!(graph, [0.0], line=:dash, color=:orange, lw=2, label="")

    if neighbour_order_to_measure_to==2
        annotation = L"E^{(1)} = E^{(2)}"
    elseif neighbour_order_to_measure_to==1
        annotation = L"E^{(0)} = E^{(1)}"
    else
        error("Invalid neighbour_order_to_measure_to value for annotation: $neighbour_order_to_measure_to")
    end

    if connectivity == "Slice-Rotation"
        annotate!(graph, [(xlims(graph)[1]+0.05, 0.0-0.003, Plots.text(annotation, 10, :black))])
    else
        annotate!(graph, [(xlims(graph)[1]+0.05, 0.0-0.0065, Plots.text(annotation, 10, :black))])
    end

    # Add E_star vertical lines to graphs if E_star is nonzero and we are slice rotation cube
    if E_star != 0.0 && connectivity == "Slice-Rotation"
        vline!(graph, [E_star], line=:dash, color=:green, lw=2, label="")
        annotate!(graph, [(E_star+0.05, ylims(graph)[1]-0.05, Plots.text(L"E^*", 8, :black))])
    end



    ### -- Save and display the graph --
    println("Saving horizontal graph...")
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram.svg")
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram.png")
    display(graph)

end

function remove_bad_rows(data::Array{Float64,2}, L::Int64)::Array{Float64,2}
    # Find rows without NaN or negative values or above solved configuration energy in the first two columns or above solved configuration energy
    non_bad_rows = .!isnan.(data[:, 1]) .& .!isnan.(data[:, 2]) .& (data[:, 1] .>= 0) .& (data[:, 2] .>= 0) .& (data[:, 1] .< -solved_configuration_energy(RubiksCube(L))) .& (data[:, 2] .< -solved_configuration_energy(RubiksCube(L)))
    # Return the data without rows containing NaN
    return data[non_bad_rows, :]
end