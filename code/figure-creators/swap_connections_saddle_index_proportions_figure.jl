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


function saddle_index_proportions_figure(N::Int64, neighbour_order_to_measure_to::Int64=1)
    connectivity = "Swap-Move"
    simulation_name = "swap_$(N)_connections"

    # E_star = -0.39015151515151514 
    # E_star = -0.376098787878788
    E_star = -0.3859651515151515
    E_on = -0.23134348484848488


    ### -- SET UP DEFAULT PARAMETERS --
    header_line = readlines(joinpath("results/neighbour_initial_and_final_energies_distribution_results/swap_connections_results/combined_data/swap_N_$(N)_combined_disorder_average_connections_L=5_inherent_disorder_E0_E1_swap_E0_E1_energy_saddle_index_densities"))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    ### --- READ IN DATA ---
    energy_saddle_index_densities_data_matrix = readdlm(joinpath("results/neighbour_initial_and_final_energies_distribution_results/swap_connections_results/combined_data/swap_N_$(N)_combined_disorder_average_connections_L=5_inherent_disorder_E0_E1_swap_E0_E1_energy_saddle_index_densities"), ',', Float64, '\n', skipstart=3)

    # Remove any rows with energy values larger than e10
    energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .< 1e10, :]
    # or with energy values >=0.0
    energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .< abs(solved_configuration_energy(L)), :]
    # or with energy values <= -solved_configuration_energy(L)
    energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[energy_saddle_index_densities_data_matrix[:,1] .> 0, :]

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

    # Z = 6*(L-1)
    Z = N

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


        # Calculate the proportion of minima
        k_saddle_proportions[i, 1] = sum(E0_saddle_indices .== 0) / length(E0_saddle_indices)
        k_saddle_proportions[i, 2] = sum(E0_saddle_indices .== 1) / length(E0_saddle_indices)
        k_saddle_proportions[i, 3] = sum(E0_saddle_indices .== 2) / length(E0_saddle_indices)
        k_saddle_proportions[i, 4] = sum(E0_saddle_indices .> 1) / length(E0_saddle_indices)
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

    scatter!(k_saddle_proportions_graph,
        -1.0 .+ (E0_bin_values./-solved_configuration_energy(cube)),
        k_saddle_proportions[:, 1],
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel="Saddle Index, "*L"K"*", Proportions",
        # title="L=$L $connectivity Cube Saddle Index Proportions",
        # titlefontsize=10,   # Title font size

        # label= neighbour_order_to_measure_to==1 ? "K=0 (Minima)" : "K=$(0+neighbour_order_to_measure_to-1)",
        label="",
        xguidefontsize=12,   # X-axis label font size
        yguidefontsize=12,   # Y-axis label font size
        margin=0mm,          # Margin around the plot
        legend=:right,
        color=alex_blue,
    )

    # Empty scatter just to put K>2 in correct legend and color position
    scatter!(k_saddle_proportions_graph, [], [], label="K≥$(2+neighbour_order_to_measure_to-1)", color=alex_green)


    if connectivity=="Slice-Rotation"

        # Add E^* vertical line
        if !isnothing(E_star)
            vline!(k_saddle_proportions_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="")
            annotate!(k_saddle_proportions_graph, [(E_star+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\epsilon^*", 12, :black))])
        end

        # Add E^on vertical line
        if !isnothing(E_on)
            vline!(k_saddle_proportions_graph, [E_on], linecolor=:red, linestyle=:dash, linewidth=2, label="")
            annotate!(k_saddle_proportions_graph, [(E_on+0.023, ylims(k_saddle_proportions_graph)[1]+0.58, Plots.text(L"\epsilon^{\rm on}", 12, :black))])
        end

    end

    # Save a-0.24010378787878783nd display the graphs
    # Remove string "_energy_saddle_index_densities" from simulation_name to get save_name
    save_name = simulation_name


    if connectivity=="Slice-Rotation"
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions.pdf")
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions.png")
    else 
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions_$(connectivity).pdf")
        savefig(k_saddle_proportions_graph, "results/neighbour_initial_and_final_energies_distribution_results/$(save_name)_k_saddle_proportions_$(connectivity).png")
    end
    
    
    display(k_saddle_proportions_graph)

end














function combined_swap_connections_saddle_index_proportions_figure(N_values::Array{Int64}; include_K1_in_K0::Bool=false)
    
    L = 5
    neighbour_order_to_measure_to::Int64=1

    addon = include_K1_in_K0 ? "_with_K1_in_K0" : ""

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)

    # Colors and symbols for different N values
    colors = [RGB(227/255, 11/255, 92/255), RGB(255/255, 165/255, 0/255), RGB(23/255,177/255,105/255), RGB(100/255, 149/255, 237/255), RGB(255/255, 105/255, 180/255)]



    # E_star = -0.39015151515151514 
    # E_star = -0.376098787878788
    E_star = -0.3859651515151515
    # E_on = -0.23134348484848488
    E_on = -0.24010378787878783


    ### --- READ IN DATA ---
    data_dict = Dict()

    for (i, N) in enumerate(N_values)
        file_path = joinpath("results/neighbour_initial_and_final_energies_distribution_results/swap_connections_results/combined_data/swap_N_$(N)_combined_disorder_average_connections_L=5_inherent_disorder_E0_E1_swap_E0_E1_energy_saddle_index_densities")
        
        cube = RubiksCube(L)
        
        energy_saddle_index_densities_data_matrix = readdlm(file_path, ',', Float64, '\n', skipstart=3)
        
        # Remove rows with energy values larger than e10 or >=0.0
        energy_saddle_index_densities_data_matrix = energy_saddle_index_densities_data_matrix[
            (energy_saddle_index_densities_data_matrix[:,1] .< 1e10) .& 
            (energy_saddle_index_densities_data_matrix[:,1] .< abs(solved_configuration_energy(N))), :]
        
        E0_bin_values = unique(energy_saddle_index_densities_data_matrix[:,1])
        
        Z = N
        normalization_factor = Z*(Z-1)^(neighbour_order_to_measure_to-1)
        neighbours_saddle_indices_shared_between = Z^(neighbour_order_to_measure_to-1)
        
        k_saddle_proportions = zeros(Float64, length(E0_bin_values), 4)
        
        for (j, E0) in pairs(E0_bin_values)
            E0_saddle_indices = [round((normalization_factor*energy_saddle_index_densities_data_matrix[k,2])/(neighbours_saddle_indices_shared_between)) 
                                 for k in 1:size(energy_saddle_index_densities_data_matrix,1) 
                                 if energy_saddle_index_densities_data_matrix[k,1] == E0]
            
            k_saddle_proportions[j, 1] = sum(E0_saddle_indices .== 0) / length(E0_saddle_indices)
            k_saddle_proportions[j, 2] = sum(E0_saddle_indices .== 1) / length(E0_saddle_indices)
            
            if include_K1_in_K0
                k_saddle_proportions[j,1] = k_saddle_proportions[j,1] + k_saddle_proportions[j,2]
            end
            
            k_saddle_proportions[j, 4] = sum(E0_saddle_indices .> 1) / length(E0_saddle_indices)
        end
        
        data_dict[N] = (E0_bin_values, k_saddle_proportions)
    end







    ### --- CREATE COMBINED GRAPH ---
    combined_plot = scatter(
        xlabel="Energy Density, "*L"\epsilon = E/|\!\!E_s|",
        ylabel=L"K=0"*" and "*L"K \geq 2"*" Proportions",
        legend=:right,
        xguidefontsize=12,
        yguidefontsize=12,
        margin=5mm,
        xlims=(-0.7, -0.05)
    )

    for (i, N) in enumerate(N_values)
        E0_bin_values, k_saddle_proportions = data_dict[N]
        
        scatter!(combined_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 1],
            label="K=0, N=$N",
            # label="",
            color=colors[i],
            marker=:circle,
            markersize=4
        )

        scatter!(combined_plot,
            -1.0 .+ (E0_bin_values./-solved_configuration_energy(RubiksCube(L))),
            k_saddle_proportions[:, 4],
            label="K≥2, N=$N",
            # label="",
            color=colors[i],
            marker=:diamond,
            markersize=4,
        )
    end


    display(combined_plot)
    savefig(combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/combined_swap_connections_saddle_index_proportions_L=5_varyng_N$(addon).pdf")
    savefig(combined_plot, "results/neighbour_initial_and_final_energies_distribution_results/combined_swap_connections_saddle_index_proportions_L=5_varyng_N$(addon).png")
end













# saddle_index_proportions_figure(100)
# saddle_index_proportions_figure(500)
# saddle_index_proportions_figure(1000)
# saddle_index_proportions_figure(2000)