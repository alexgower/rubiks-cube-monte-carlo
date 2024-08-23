using DelimitedFiles
using LaTeXStrings, Plots; 

using StatsBase
using CSV
# using DataFrames

using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

using Images

include("../core/rubiks_cube.jl")

# Make nice function renamings
E1_E2_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation"; bin_diagonal_graph=false) = energy_connectivity_histogram_figure(simulation_name, connectivity; neighbour_order_to_measure_to=2, bin_diagonal_graph=bin_diagonal_graph)
E0_E1_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation"; bin_diagonal_graph=false) = energy_connectivity_histogram_figure(simulation_name, connectivity; neighbour_order_to_measure_to=1, bin_diagonal_graph=bin_diagonal_graph)

function annotatewithbox!(
    fig::Plots.Plot,
    text::Plots.PlotText,
    x::Real, y::Real,
    Δx::Real, Δy::Real = Δx; kwargs...)
    
    box = Plots.Shape(:rect)
    Plots.scale!(box, Δx, Δy)
    Plots.translate!(box, x, y)
    
    Plots.plot!(fig, box, c=:white, linecolor=:white, label=false; kwargs...)
    Plots.annotate!(fig, x+Δx/30, y, text)
    
    fig
end

# Main function
function energy_connectivity_histogram_figure(simulation_name::String, connectivity::String="Slice-Rotation"; neighbour_order_to_measure_to::Int64=1, bin_diagonal_graph::Bool=false, bin_horizontal_graph::Bool=false, prediction::Bool=false)

    # E_star = -0.39015151515151514
    E_star = -0.376098787878788

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
        # Print the slice sum for each E_0 slice
        println("E_0 = $(bin_edges_x[i]): $(slice_sum)")

        # Also print the proportion of connections between +1 and +10 avoe this E_0
        mask = (energy_connections_data_matrix[:,1] .== bin_edges_x[i]) .& (energy_connections_data_matrix[:,2] .> bin_edges_x[i] .+ 1) .& (energy_connections_data_matrix[:,2] .<= bin_edges_x[i] .+ 10)
        proportion = sum(mask) / slice_sum
        println("Proportion of connections between +1 and +10 above this E_0: $(proportion)")

    end

    E0_values = energy_connections_data_matrix[:,1]./-solved_configuration_energy(cube)
    E1_values = energy_connections_data_matrix[:,2]./-solved_configuration_energy(cube)

    if bin_diagonal_graph==false

    if neighbour_order_to_measure_to == 2
        # diagonal_ylabel = "Second-Neighbour Configuration Energy Density, "*L"E^{(2)}/E_s"
        # diagonal_xlabel = "Neighbour Configuration Energy Density, "*L"E^{(1)}/E_s"
        diagonal_ylabel = "Second-Neighbour Energy Density, "*L"\epsilon^{(2)}"
        diagonal_xlabel = "Neighbour Energy Density, "*L"\epsilon^{(1)}"
    elseif neighbour_order_to_measure_to == 1
        # diagonal_ylabel = "Neighbour Configuration Energy Density, "*L"E^{(1)}/E_s"
        # diagonal_xlabel = "Initial Configuration Energy Density, "*L"E^{(0)}/E_s"
        diagonal_ylabel = "Neighbour Energy Density, "*L"\epsilon^{(1)}"
        diagonal_xlabel = "Energy Density, "*L"\epsilon^{(0)}"
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
        # title="$connectivity Cube", 
        xlims=(minimum(E0_values), maximum(E0_values)), 
        ylims=(minimum(E1_values), maximum(E1_values)),
        colorbar_title="",
        # titlefontsize=10,   # Title font size
        xguidefontsize=12,   # X-axis label font size
        yguidefontsize=12,   # Y-axis label font size
        margin=5mm,          # Margin around the plot
        xticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
        yticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
        # background_color=:white,
        # grid=false
    )
    annotate!(1.115 * maximum(E0_values), 0.5 * (minimum(E1_values) + maximum(E1_values)), text("Sampled Frequency", 10, :center, :center, rotation=90))

    # Add E_0 = E1 lines to graph
    if neighbour_order_to_measure_to==2
        # annotation = L"E⁽²⁾=E⁽¹⁾"
        annotation = L"\epsilon^{(2)} = \epsilon^{(1)}"
    elseif neighbour_order_to_measure_to==1
        # annotation = L"E⁽¹⁾=E⁽⁰⁾"
        annotation = L"\epsilon^{(1)} = \epsilon^{(0)}"
    else
        error("Invalid neighbour_order_to_measure_to value for annotation: $neighbour_order_to_measure_to")
    end

    annotate!(graph, [(xlims(graph)[1]+0.1, ylims(graph)[1]+0.03, Plots.text(annotation, 12, :black))])

    plot!(graph, [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], [min_value/-solved_configuration_energy(cube), max_value/-solved_configuration_energy(cube)], line=:dash, color=:orange, lw=2, label="")

    # Add E_star vertical lines to graphs if we are slice rotation cube
    if connectivity == "Slice-Rotation"
        E_star_plot = 1+E_star
        vline!(graph, [E_star_plot], line=:dash, color=:green, lw=2, label="")
        annotate!(graph, [(E_star_plot+0.02, ylims(graph)[1]+0.05, Plots.text(L"\epsilon^*", 12, :black))])

    end

        # Add title as annotated text in top right corner
        if connectivity=="Slice-Rotation"
            annotate!(graph, [(xlims(graph)[1]+0.12, ylims(graph)[2]-0.25, Plots.text("$(connectivity) Cube", 10, :black))])
        else
            annotate!(graph, [(xlims(graph)[2]-0.1, ylims(graph)[2]-0.25, Plots.text("$(connectivity) Cube", 10, :black))])
        end

    # ### --- PLOT SPIRAL IMAGES ON GRAPH ---
    # if connectivity == "Slice-Rotation" && neighbour_order_to_measure_to == 1
    #     ## -- SLICE SPIRAL --
    #     slice_spiral_image = load("results/final_paper_results/spiral.png")

    #     # Determine the desired width and height on the graph
    #     # Here you set one dimension, and calculate the other to preserve the aspect ratio
    #     desired_width = 0.4
    #     aspect_ratio = size(slice_spiral_image, 2) / size(slice_spiral_image, 1) # width / height
    #     desired_height = (desired_width / aspect_ratio)

    #     # Determine the location on the graph where you want the image's bottom-left corner
    #     x_location = 0.13
    #     y_location = 0.065

    #     # Plot the image with the specified dimensions and location
    #     plot!(graph, reverse(slice_spiral_image; dims=1), yflip=false, inset=bbox(x_location,y_location,desired_width-0.1,desired_height), subplot=2, aspect_ratio=:auto, axis=false, grid=false, framestyle=:box, legend=false, ticks=nothing, border=:none, plot_bgcolor=:transparent)


    #     ## -- SWAP SPIRAL --
    #     swap_spiral_image = load("results/final_paper_results/swap-spiral.png")

    #     # Determine the desired width and height on the graph
    #     # Here you set one dimension, and calculate the other to preserve the aspect ratio
    #     desired_width = 0.4
    #     aspect_ratio = size(swap_spiral_image, 2) / size(swap_spiral_image, 1) # width / height
    #     desired_height = (desired_width / aspect_ratio)

    #     # Determine the location on the graph where you want the image's bottom-left corner
    #     x_location = 0.73
    #     y_location = 0.5

    #     # Plot the image with the specified dimensions and location
    #     plot!(graph, reverse(swap_spiral_image; dims=1), yflip=false, inset=bbox(x_location,y_location,desired_width-0.1,desired_height), subplot=2, aspect_ratio=:auto, axis=false, grid=false, framestyle=:box, legend=false, ticks=nothing, border=:none, plot_bgcolor=:transparent)

    # end

    ### -- Save and display the graphs --
    println("Saving diagonal graph...")
    if connectivity=="Slice-Rotation"
        extra = prediction ? "_prediction" : ""

        savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal_raw$(extra).pdf")
        savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal_raw$(extra).png")
        display(graph)
    else
        savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal.pdf")
        savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram_diagonal.png")
        display(graph)
    end
    

    end


























    ## -- HORIZONTAL GRAPH --
    if bin_horizontal_graph
        return
    end

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
        # horizontal_ylabel = "Second-Neighbour Energy Density Difference, "*L"(E^{(2)}\!\!\!-\!E^{(1)})/\,|\!\!E_s\,|"
        horizontal_ylabel = "Second-Neighbour Energy Density Difference, "*L"\epsilon^{(2)}-\epsilon^{(1)}"        
        # horizontal_xlabel = "Neighbour Ernergy Density, "*L"E^{(1)}\!\!/\,|\!\!E_s\,|"
        horizontal_xlabel = "Neighbour Energy Density, "*L"\epsilon^{(1)}"
    elseif neighbour_order_to_measure_to == 1
        # horizontal_ylabel = "Neighbour Energy Density Difference, "*L"(E^{(1)}\!\!\!-\!E^{(0)})/\,|\!\!E_s\,|"
        horizontal_ylabel = "Neighbour Energy Density Difference, "*L"\epsilon^{(1)}-\epsilon^{(0)}"
        # horizontal_xlabel = "Energy Density, "*L"E^{(0)}\!\!/\,|\!\!E_s\,|"
         horizontal_xlabel = "Energy Density, "*L"\epsilon^{(0)}"
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
        colorbar_title="",
        xguidefontsize=8,   # X-axis label font size
        yguidefontsize=8,   # Y-axis label font size
        margin=5mm,          # Margin around the plot
        # title="$connectivity Cube",
        # titlefontsize=10,   # Title font size
        xticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
    )
    annotate!(1.115 * maximum(E0_values), 0.5 * (minimum(E_difference_values) + maximum(E_difference_values)), text("Sampled Frequency", 8, :center, :center, rotation=90))

    # Add title as annotated text in top right corner
    if connectivity=="Slice-Rotation"
        annotate!(graph, [(xlims(graph)[2]-0.1, ylims(graph)[2]-0.012, Plots.text("$(connectivity) Cube", 10, :black))])
    else
        annotate!(graph, [(xlims(graph)[2]-0.1, ylims(graph)[2]-0.012, Plots.text("$(connectivity) Cube", 10, :black))])
    end

    # Add E_0 = E1 lines to graph
    hline!(graph, [0.0], line=:dash, color=:orange, lw=2, label="")

    if neighbour_order_to_measure_to==2
        annotation = L"\epsilon^{(2)} = \epsilon^{(1)}"
    elseif neighbour_order_to_measure_to==1
        annotation = L"\epsilon^{(1)} = \epsilon^{(0)}" 
    else
        error("Invalid neighbour_order_to_measure_to value for annotation: $neighbour_order_to_measure_to")
    end

    if connectivity == "Slice-Rotation"
        if neighbour_order_to_measure_to==1
            annotate!(graph, [(xlims(graph)[1]+0.05, 0.0-0.005, Plots.text(annotation, 10, :black))])
        else
            xpos = xlims(graph)[1]+0.1
            ypos = 0.0-0.0065
            rectangle_width = 0.06
            rectangle_height = 0.004
            annotatewithbox!(graph, Plots.text(annotation, 10, color), xpos, ypos, rectangle_width, rectangle_height)

        end
    else
        annotate!(graph, [(xlims(graph)[1]+0.05, 0.0-0.010, Plots.text(annotation, 10, :black))])
        # annotate!(graph, [(xlims(graph)[1]+0.05, 0.0-0.0115, Plots.text(annotation, 10, :black))])
    end


    # Add E_star vertical lines to graphs if we are slice rotation cube
    if connectivity == "Slice-Rotation"
        E_star_plot = 1+E_star
        vline!(graph, [E_star_plot], line=:dash, color=:green, lw=2, label="")
        annotate!(graph, [(E_star_plot+0.02, ylims(graph)[1]+0.005, Plots.text(L"\epsilon^*", 10, :black))])

        if prediction
            x = range(minimum(E0_values), stop=maximum(E0_values), length=100) 
            energy_densities = x .- 1
            # Plot the function (+|\epsilon|*(2*4*L) - (1/6)*(2*4*11))/|solved_configuration_energy(cube)|
            y = (abs.(energy_densities).*(2*4*L) .- (1/6)*(2*4*11)) ./ abs(solved_configuration_energy(cube))
            plot!(graph, x, y, line=:dash, color=:red, lw=2, label="")

            y = (abs.(energy_densities).*(4*L) .- (1/6)*(4*11)) ./ abs(solved_configuration_energy(cube))
            plot!(graph, x, y, line=:dash, color=:pink, lw=2, label="")

        end
    end



    ### -- Save and display the graph --
    println("Saving horizontal graph...")
    extra = prediction ? "_prediction" : ""
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram$(extra).pdf")
    savefig(graph, "results/final_paper_results/$(simulation_name)_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_histogram$(extra).png")
    display(graph)

end

function remove_bad_rows(data::Array{Float64,2}, L::Int64)::Array{Float64,2}
    # Find rows without NaN or negative values or above solved configuration energy in the first two columns or above solved configuration energy
    non_bad_rows = .!isnan.(data[:, 1]) .& .!isnan.(data[:, 2]) .& (data[:, 1] .>= 0) .& (data[:, 2] .>= 0) .& (data[:, 1] .< -solved_configuration_energy(RubiksCube(L))) .& (data[:, 2] .< -solved_configuration_energy(RubiksCube(L)))
    # Return the data without rows containing NaN
    return data[non_bad_rows, :]
end