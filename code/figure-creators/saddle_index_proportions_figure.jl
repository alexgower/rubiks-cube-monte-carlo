using Pkg
Pkg.activate("/home/apg59/rubiks-cube-monte-carlo")

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


function saddle_index_proportions_figure(connectivity="Slice-Rotation", neighbour_order_to_measure_to::Int64=1)
    simulation_name = "L=11_combined"

    simulation_name = simulation_name*"_energy_saddle_index_densities_slice"

    # E_star = -0.39015151515151514 
    # E_star = -0.376098787878788
    E_star = -0.3859651515151515
    # E_on = -0.23134348484848488
    E_on = -0.24010378787878783


    ### -- SET UP DEFAULT PARAMETERS --
    header_line = readlines(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_disorder_average_connections_L=11_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    ### --- READ IN DATA ---
    energy_saddle_index_densities_data_matrix = readdlm(joinpath("results/neighbour_initial_and_final_energies_distribution_results/E1_E0_results/combined_data/combined_disorder_average_connections_L=11_inherent_disorder_E0_E1_slice_E0_E1_energy_saddle_index_densities"), ',', Float64, '\n', skipstart=3)

    # Remove any rows with energy values larger than e10
    energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[abs.(energy_saddle_index_densities_data_matrix[:,1]) .< 1e10, :]
    # or with energy values >=0.0
    energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .< abs(solved_configuration_energy(L)), :]
    # or with energy values < 0
    energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .> 0, :]


    ### --- COLOURS ---
    Plots.default(dpi = 600)

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
    k_saddle_errors = zeros(Float64, length(E0_bin_values), 4)  # New error array

    Z = 6*(L-1)
    normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
    neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)

    K_60 = 0
    K_0 = 0

    total = 0
    for (i, E0) in pairs(E0_bin_values)
        # Get the saddle indices for the current E0 slice
        E0_saddle_indices = [round((normalization_factor*energy_saddle_index_densities_data_matrix[j,2])/(neighbours_saddle_indices_shared_between)) for j in 1:size(energy_saddle_index_densities_data_matrix,1) if energy_saddle_index_densities_data_matrix[j,1] == E0] 
        
        K_60 += sum(E0_saddle_indices .== 60)
        K_0 += sum(E0_saddle_indices .== 0)
        total += length(E0_saddle_indices)

        N = length(E0_saddle_indices)

        # Calculate the proportion of minima
        k_saddle_proportions[i, 1] = sum(E0_saddle_indices .== 0) / N
        k_saddle_proportions[i, 2] = sum(E0_saddle_indices .== 1) / N
        k_saddle_proportions[i, 3] = sum(E0_saddle_indices .== 2) / N
        k_saddle_proportions[i, 4] = sum(E0_saddle_indices .> 1) / N
        
        # Calculate standard errors for proportions
        for j in 1:4
            p = k_saddle_proportions[i, j]
            k_saddle_errors[i, j] = sqrt(p * (1 - p) / N)
        end
    end

    println("K=60 Maxima Saddles: ", K_60)
    println("K=0 Minima Saddles: ", K_0)
    println("Total Saddles: ", total)

    # Plot the k saddle proportions against E0
    # Empty scatter just to put K=0 in correct legend and color position
    k_saddle_proportions_graph = scatter([], [],         label= neighbour_order_to_measure_to==1 ? "K=0 (Minima)" : "K=$(0+neighbour_order_to_measure_to-1)", color=alex_blue)

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
    k_saddle_proportions[:, 4],
    yerror=k_saddle_errors[:, 4],
    label="",
    color=alex_green,
    )

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 2],
        yerror=k_saddle_errors[:, 2],
        label="K=$(1+neighbour_order_to_measure_to-1)",
        color=alex_pink,
    )

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 3],
        yerror=k_saddle_errors[:, 3],
        label="K=$(2+neighbour_order_to_measure_to-1)",
        color=alex_orange,
    )

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 1],
        yerror=k_saddle_errors[:, 1],
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Saddle Index, "*L"K"*", Proportions, "*L"\bar p_{K}(\epsilon)",
        # title="L=$L $connectivity Cube Saddle Index Proportions",
        # titlefontsize=10,   # Title font size

        # label= neighbour_order_to_measure_to==1 ? "K=0 (Minima)" : "K=$(0+neighbour_order_to_measure_to-1)",
        label="",
        xguidefontsize=12,   # X-axis label font size
        yguidefontsize=12,   # Y-axis label font size
        margin=1mm,          # Margin around the plot
        legend=:left,
        color=alex_blue,
    )

    # Empty scatter just to put K>2 in correct legend and color position
    scatter!(k_saddle_proportions_graph, [], [], label="K≥$(2+neighbour_order_to_measure_to-1)", color=alex_green)


    if connectivity=="Slice-Rotation"

        # Add E^* vertical line
        if !isnothing(E_star)
            vline!(k_saddle_proportions_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="")
            annotate!(k_saddle_proportions_graph, [(E_star+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\bar\epsilon^*", 12, :black))])
        end

        # Add E^on vertical line
        if !isnothing(E_on)
            vline!(k_saddle_proportions_graph, [E_on], linecolor=alex_red, linestyle=:dash, linewidth=2, label="")
            annotate!(k_saddle_proportions_graph, [(E_on+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\bar\epsilon^{\rm on}", 12, :black))])
        end

    end

    # Save and display the graphs
    # Remove string "_energy_saddle_index_densities" from simulation_name to get save_name
    save_name = replace(simulation_name, "_energy_saddle_index_densities" => "")


    if connectivity=="Slice-Rotation"
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions.pdf")
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions.png")
    else 
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions_$(connectivity).pdf")
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions_$(connectivity).png")
    end
    
    
    display(k_saddle_proportions_graph)

end


saddle_index_proportions_figure()