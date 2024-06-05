using DelimitedFiles
using Plots

# using StatsBase
using LaTeXStrings
using CSV
# using DataFrames

# using Colors
# using ColorTypes, ColorSchemes
using Plots.PlotMeasures

include("../core/rubiks_cube.jl")


function saddle_index_proportions_figure(simulation_name::String; neighbour_order_to_measure_to::Int64=1)

    # E_star = -0.39015151515151514 
    E_star = -0.376098787878788


    ### -- SET UP DEFAULT PARAMETERS --
    header_line = readlines(joinpath("results/final_paper_results",simulation_name))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    ### --- READ IN DATA ---
    energy_saddle_index_densities_data_matrix = readdlm(joinpath("results/final_paper_results",simulation_name), ',', Float64, '\n', skipstart=3)

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)


    ### --- CREATE GRAPH ---
    E0_bin_values = unique(energy_saddle_index_densities_data_matrix[:,1])

    
    # Now gather the proportion of minima (saddle index K=0), K=1, K=2  and K>2 for each energy value
    # and store in a new matrix against energy
    k_saddle_proportions = zeros(Float64, length(E0_bin_values), 4)

    Z = 6*(L-1)
    normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
    neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)

    for (i, E0) in pairs(E0_bin_values)
        # Get the saddle indices for the current E0 slice
        E0_saddle_indices = [round((normalization_factor*energy_saddle_index_densities_data_matrix[j,2])/(neighbours_saddle_indices_shared_between)) for j in 1:size(energy_saddle_index_densities_data_matrix,1) if energy_saddle_index_densities_data_matrix[j,1] == E0] 
        # Calculate the proportion of minima
        k_saddle_proportions[i, 1] = sum(E0_saddle_indices .== 0) / length(E0_saddle_indices)
        k_saddle_proportions[i, 2] = sum(E0_saddle_indices .== 1) / length(E0_saddle_indices)
        k_saddle_proportions[i, 3] = sum(E0_saddle_indices .== 2) / length(E0_saddle_indices)
        k_saddle_proportions[i, 4] = sum(E0_saddle_indices .> 1) / length(E0_saddle_indices)
    end

    # Plot the k saddle proportions against E0
    k_saddle_proportions_graph = scatter(
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 1],
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Saddle Index, "*L"K"*", Proportions",
        # title="L=$L $connectivity Cube Saddle Index Proportions",
        # titlefontsize=10,   # Title font size

        label= neighbour_order_to_measure_to==1 ? "K=0 (Minima)" : "K=$(0+neighbour_order_to_measure_to-1)",
        xguidefontsize=12,   # X-axis label font size
        yguidefontsize=12,   # Y-axis label font size
        margin=0mm,          # Margin around the plot
        legend=:left,
        color=alex_blue,
    )

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
    k_saddle_proportions[:, 4],
    label="",
    color=alex_green,
    )

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 2],
        label="K=$(1+neighbour_order_to_measure_to-1)",
        color=alex_pink,
    )

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 3],
        label="K=$(2+neighbour_order_to_measure_to-1)",
        color=alex_orange,
    )

    # Empty scatter just to put K>2 in correct legend and color position
    scatter!(k_saddle_proportions_graph, [], [], label="K≥$(2+neighbour_order_to_measure_to-1)", color=alex_green)



    
    # Add E^* vertical line
    if !isnothing(E_star)
        vline!(k_saddle_proportions_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="")
        annotate!(k_saddle_proportions_graph, [(E_star+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\epsilon^*", 12, :black))])
    end

    # Save and display the graphs
    # Remove string "_energy_saddle_index_densities" from simulation_name to get save_name
    save_name = replace(simulation_name, "_energy_saddle_index_densities" => "")

    println("results/final_paper_results/$(save_name)_k_saddle_proportions.svg")

    savefig(k_saddle_proportions_graph, "results/final_paper_results/$(save_name)_k_saddle_proportions.svg")
    savefig(k_saddle_proportions_graph, "results/final_paper_results/$(save_name)_k_saddle_proportions.png")
    display(k_saddle_proportions_graph)
end
