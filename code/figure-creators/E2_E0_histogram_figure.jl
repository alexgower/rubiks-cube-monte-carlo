using DelimitedFiles
using LaTeXStrings, Plots; 

using StatsBase
using CSV
# using DataFrames

using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

# using Images

include("../core/rubiks_cube.jl")


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
function E2_E0_histogram_figure()

    connectivity = "Slice-Rotation"

    ### --- COLOURS ---
    Plots.default(dpi = 600)

    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)
    alex_grey = RGB(113/255, 121/255, 126/255)

    alex_alt_blue = RGB(4/255, 57/255, 94/255)

    # E_star = -0.39015151515151514
    # E_star = -0.376098787878788
    E_star = -0.3859651515151515
    # E_on = -0.23134348484848488
    E_on = -0.24010378787878783


    ### --- READ IN THE DATA ---
    # filename = "results/neighbour_initial_and_final_energies_distribution_results/E2_E0_results",simulation_name*"_E$(neighbour_order_to_measure_to-1)_E$(neighbour_order_to_measure_to)_energy_connections"
    filename = "results/neighbour_initial_and_final_energies_distribution_results/E2_E0_results/combined_data/combined_disorder_average_connections_L=11_inherent_disorder_E0_E1_E2_slice_energy_connections"

    header_line = readlines(joinpath(filename))[1]
    energy_connections_data_matrix = readdlm(joinpath(filename), ',', Float64, '\n', skipstart=2)


    ### --- SET UP DEFAULT PARAMETERS ---

    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    # Remove NaNs and overflow # TODO how exist?
    energy_connections_data_matrix = remove_bad_rows(energy_connections_data_matrix, L)





    # ### --- HISTOGRAM GENERAL CALCULATIONS ---

    # # Determine the bin edges based on the data
    # bin_edges_x = 0:1:-solved_configuration_energy(cube)
    # bin_edges_y = 0:1:-solved_configuration_energy(cube)
    # edges = (bin_edges_x, bin_edges_y)

    # # Compute the 2D histogram
    # hist_2d = fit(Histogram, (energy_connections_data_matrix[:,1], energy_connections_data_matrix[:,2]), edges)


    # ### -- Plots. jl 2D HISTOGRAM --

    # min_value = minimum([minimum(energy_connections_data_matrix[:,1]), minimum(energy_connections_data_matrix[:,2])])
    # max_value = maximum([maximum(energy_connections_data_matrix[:,1]), maximum(energy_connections_data_matrix[:,2])])


    # # Initialize an array to hold the biases
    # biases = zeros(size(energy_connections_data_matrix, 1))

    # # NOTE THAT BIASES BELOW DO NOT MEAN FREQUENCIES IN THE HISTOGRAM BUT RATHER HOW WE BIAS THESE FREQUENCIES
    # # Calculate the sum of counts in each E_0 slice
    # for i in 1:size(hist_2d.weights, 1)
    #     slice_sum = sum(hist_2d.weights[i, :])
    #     if slice_sum != 0  # Avoid division by zero
    #         # Assign the reciprocal of the slice sum to the weights array for each data point in the current E_0 slice
    #         mask = energy_connections_data_matrix[:,1] .== bin_edges_x[i]
    #         biases[mask] .= 1.0 / slice_sum
    #     end
    #     # Print the slice sum for each E_0 slice
    #     println("E_0 = $(bin_edges_x[i]): $(slice_sum)")

    #     # Also print the proportion of connections between +1 and +10 avoe this E_0
    #     # mask = (energy_connections_data_matrix[:,1] .== bin_edges_x[i]) .& (energy_connections_data_matrix[:,2] .> bin_edges_x[i] .+ 1) .& (energy_connections_data_matrix[:,2] .<= bin_edges_x[i] .+ 10)
    #     # proportion = sum(mask) / slice_sum
    #     # println("Proportion of connections between +1 and +10 above this E_0: $(proportion)")

    # end






    ## -- HORIZONTAL GRAPH --

    # Determine the bin edges based on the data
    bin_edges_x = 0:1:-solved_configuration_energy(cube)
    bin_edges_y = solved_configuration_energy(cube):1:-solved_configuration_energy(cube)
    edges = (bin_edges_x, bin_edges_y)

    # Compute the 2D histogram
    hist_2d = fit(Histogram, (energy_connections_data_matrix[:,1], energy_connections_data_matrix[:,3].-energy_connections_data_matrix[:,1]), edges)


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
    E_difference_values = (energy_connections_data_matrix[:,3] .- energy_connections_data_matrix[:,1])./-solved_configuration_energy(cube)

    horizontal_ylabel = L"\epsilon^{(2)}-\epsilon^{(0)}"        
    horizontal_xlabel = L"\epsilon^{(0)}"


    graph = histogram2d(
        E0_values, 
        E_difference_values, 
        # color=:bluesreds,
        # color=cgrad([alex_blue, :white], [0.0, 0.15, 0.3]),
        color = connectivity == "Slice-Rotation" ? cgrad([:darkblue, :white], [0.0, 0.15, 0.3]) : cgrad([:darkblue, :white], [0.0, 0.05, 0.1]),        show_empty_bins=false,
        # normalize=:pdf, 
        bins=((bin_edges_x./-solved_configuration_energy(cube)), bin_edges_y./-solved_configuration_energy(cube)), 
        weights=biases, 
        ylabel=horizontal_ylabel, 
        xlabel=horizontal_xlabel,
        xlims=(minimum(E0_values)+0.02, maximum(E0_values)), # Custom cleaned up
        ylims=(minimum(E_difference_values), maximum(E_difference_values)),
        colorbar_title="",
        xguidefontsize=12,   # X-axis label font size
        yguidefontsize=12,   # Y-axis label font size
        margin=5mm,          # Margin around the plot
        # title="$connectivity Cube",
        # titlefontsize=10,   # Title font size
        xticks=(0.4:0.1:0.8, -1.0 .+ 0.4:0.1:0.8),
    )
    # annotate!(1.115 * maximum(E0_values), 0.5 * (minimum(E_difference_values) + maximum(E_difference_values)), text("Sampled Frequency", 8, :center, :center, rotation=90))

    # Add title as annotated text in top right corner
    annotate!(graph, [(xlims(graph)[2]-0.35, ylims(graph)[2]-0.012, Plots.text("Slice-Rotation Cube", 10, :black))])

    # Add E_0 = E2 lines to graph
    hline!(graph, [0.0], line=:dash, color=:orange, lw=2, label="")

    annotation = L"\epsilon^{(2)} = \epsilon^{(0)}"

    annotate!(graph, [(xlims(graph)[1]+0.05, 0.0-0.005, Plots.text(annotation, 10, :black))])


    # Add E_star vertical lines to graphs if we are slice rotation cube

    E_star_plot = 1+E_star
    vline!(graph, [E_star_plot], line=:dash, color=:green, lw=2, label="")
    annotate!(graph, [(E_star_plot+0.03, ylims(graph)[1]+0.005, Plots.text(L"\bar\epsilon^*", 10, :black))])

    E_on_plot = 1+E_on
    vline!(graph, [E_on_plot], line=:dash, color=alex_red, lw=2, label="")
    annotate!(graph, [(E_on_plot-0.03, ylims(graph)[1]+0.005, Plots.text(L"\bar\epsilon^{\rm on}", 10, :black))])




    ### -- Save and display the graph --
    println("Saving horizontal graph...")
    savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/E2_E0_histogram_L=11.pdf")
    savefig(graph, "results/neighbour_initial_and_final_energies_distribution_results/E2_E0_histogram_L=11.png")
    display(graph)

end

function remove_bad_rows(data::Array{Float64,2}, L::Int64)::Array{Float64,2}
    # Find rows without NaN or negative values or above solved configuration energy in the first two columns or above solved configuration energy
    non_bad_rows = .!isnan.(data[:, 1]) .& .!isnan.(data[:, 2]) .& (data[:, 1] .>= 0) .& (data[:, 2] .>= 0) .& (data[:, 1] .< -solved_configuration_energy(RubiksCube(L))) .& (data[:, 2] .< -solved_configuration_energy(RubiksCube(L)))
    # Return the data without rows containing NaN
    return data[non_bad_rows, :]
end


E2_E0_histogram_figure()