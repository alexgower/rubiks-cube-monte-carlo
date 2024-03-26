using LaTeXStrings
using DelimitedFiles
using Plots
using StatsBase
using Plots.PlotMeasures
using Colors

include("../core/rubiks_cube.jl")

function saddle_index_density_figure(simulation_name::String)


    ### --- SET UP DEFAULT PARAMETERS ---
    header_line = readlines(joinpath("results/final_paper_results",simulation_name*"_energy_saddle_index_densities_slice"))[1]
    match_obj = match(r"L=(\d+)", header_line)
    L = parse(Int, match_obj.captures[1])
    cube = RubiksCube(L)

    E_star = -0.39015151515151514

    ### --- READ IN DATA ---
    energy_saddle_index_densities_data_matrix_slice = readdlm(joinpath("results/final_paper_results",simulation_name*"_energy_saddle_index_densities_slice"), ',', Float64, '\n', skipstart=3)    
    energy_saddle_index_densities_data_matrix_swap = readdlm(joinpath("results/final_paper_results",simulation_name*"_energy_saddle_index_densities_swap"), ',', Float64, '\n', skipstart=3)

    ### --- COLOURS ---
    alex_red = RGB(227/255, 11/255, 92/255)
    alex_pink = RGB(255/255, 105/255, 180/255)
    alex_orange = RGB(255/255, 165/255, 0/255)
    alex_green = RGB(23/255,177/255,105/255) # RGB(159/255, 226/255, 191/255)
    alex_blue = RGB(100/255, 149/255, 237/255)

    ### --- CALCULATIONS ---

    ## -- SLICE CALCULATION -- 

    # Get the unique E0 values
    E0_bin_values_slice = unique(energy_saddle_index_densities_data_matrix_slice[:,1])

    # For each E0 value, find all saddle index densities and average them
    average_saddle_index_densities_slice_including_K_1_slice = zeros(Float64, length(E0_bin_values_slice))
    average_saddle_index_densities_slice_excluding_K_1_slice = zeros(Float64, length(E0_bin_values_slice))
    for (i, E0) in pairs(E0_bin_values_slice)
        # Get the saddle index densities for the current E0 slice
        saddle_index_densities_including_K_1_slice = [energy_saddle_index_densities_data_matrix_slice[j,2] for j in 1:size(energy_saddle_index_densities_data_matrix_slice,1) if energy_saddle_index_densities_data_matrix_slice[j,1] == E0] 
        average_saddle_index_densities_slice_including_K_1_slice[i] = mean(saddle_index_densities_including_K_1_slice)

        # Now do the same calculation but if we exclude K=1 saddles
        # i.e. Remove all values whose * 6(L-1) rounds to 1.0 (i.e. K=1 saddle index)
        saddle_index_densities_excluding_K_1_slice = saddle_index_densities_including_K_1_slice[round.(saddle_index_densities_including_K_1_slice*6*(L-1)) .!= 1.0]
        average_saddle_index_densities_slice_excluding_K_1_slice[i] = mean(saddle_index_densities_excluding_K_1_slice)

    end


    ## -- SWAP CALCULATION --

    # Get the unique E0 values
    E0_bin_values_swap = unique(energy_saddle_index_densities_data_matrix_swap[:,1])

    # For each E0 value, find all saddle index densities and average them
    average_saddle_index_densities_swap = zeros(Float64, length(E0_bin_values_swap))
    for (i, E0) in pairs(E0_bin_values_swap)
        # Get the saddle index densities for the current E0 swap
        saddle_index_densities_swap = [energy_saddle_index_densities_data_matrix_swap[j,2] for j in 1:size(energy_saddle_index_densities_data_matrix_swap,1) if energy_saddle_index_densities_data_matrix_swap[j,1] == E0] 
        average_saddle_index_densities_swap[i] = mean(saddle_index_densities_swap)
    end


    ### --- PLOTTING ---

    ## -- SWAP PLOTTING --

    # Plot the average saddle index densities against E0
    # But plot all values where the average saddle index = 0.0 in green
    zero_saddle_index_densities_swap_indices = [i for (i, saddle_index_density) in pairs(average_saddle_index_densities_swap) if saddle_index_density == 0.0]
    non_zero_saddle_index_densities_indices_swap = [i for i in 1:length(E0_bin_values_swap) if average_saddle_index_densities_swap[i] != 0.0]
    println(zero_saddle_index_densities_swap_indices)

    average_saddle_index_densities_graph = plot(
        legend=:topleft,
        xguidefontsize=12,  # X-axis label font size
        yguidefontsize=12,  # Y-axis label font size
        margin=1mm,         # Margin around the plot
        xlabel="Energy Density, "*L"E/|\!\!E_s|", 
        ylabel="Average Saddle Index Density, "*L"\langle k \rangle", 
    )


    scatter!(average_saddle_index_densities_graph,
    -1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]/-solved_configuration_energy(cube)), 
    average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap], 
    label="Swap-Move Cube",
    color=alex_blue,
    )

    scatter!(average_saddle_index_densities_graph,
    -1.0 .+ (E0_bin_values_swap[zero_saddle_index_densities_swap_indices]/-solved_configuration_energy(cube)), 
    zeros(Float64, length(zero_saddle_index_densities_swap_indices)), 
    label="",
    color=alex_green,
    )

    ## -- SLICE PLOTTING --

    # Plot the average saddle index densities against E0
    # But plot all values where the average saddle index = 0.0 in red
    # And plot all values where the average saddle index = 0.0 if you exclude K=1 saddles in orange
    zero_saddle_index_densities_including_K_1_indices = [i for (i, saddle_index_density_including_K_1) in pairs(average_saddle_index_densities_slice_including_K_1_slice) if saddle_index_density_including_K_1 == 0.0]
    zero_saddle_index_densities_excluding_K_1_indices = [i for (i, saddle_index_density_excluding_K_1) in pairs(average_saddle_index_densities_slice_excluding_K_1_slice) if saddle_index_density_excluding_K_1 == 0.0 && average_saddle_index_densities_slice_including_K_1_slice[i] != 0.0]
    non_zero_saddle_index_densities_indices_slice = [i for i in 1:length(E0_bin_values_slice) if average_saddle_index_densities_slice_including_K_1_slice[i] != 0.0 && average_saddle_index_densities_slice_excluding_K_1_slice[i] != 0.0]

    println("Slice zero saddle index densities including K=1 indices: ", zero_saddle_index_densities_including_K_1_indices)
    println("Slice zero saddle index densities excluding K=1 indices: ", zero_saddle_index_densities_excluding_K_1_indices)

    scatter!(average_saddle_index_densities_graph,
        -1.0 .+ (E0_bin_values_slice[zero_saddle_index_densities_including_K_1_indices]/-solved_configuration_energy(cube)), 
        zeros(Float64, length(zero_saddle_index_densities_including_K_1_indices)), 
        label="",
        color=alex_red,
    )

    scatter!(average_saddle_index_densities_graph,
        -1.0 .+ (E0_bin_values_slice[zero_saddle_index_densities_excluding_K_1_indices]/-solved_configuration_energy(cube)), 
        zeros(Float64, length(zero_saddle_index_densities_excluding_K_1_indices)), 
        label="",
        color=alex_pink,
    )

    scatter!(average_saddle_index_densities_graph,
        -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube)), 
        average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice], 
        label="Slice-Rotation Cube",
        color=alex_orange,
    )

    # Add E^* vertical line
    vline!(average_saddle_index_densities_graph, [E_star], linecolor=:green, linestyle=:dash, linewidth=2, label="")
    annotate!(average_saddle_index_densities_graph, [(E_star+0.025, ylims(average_saddle_index_densities_graph)[1]+0.66, Plots.text(L"E^*", 12, :black))])
    

    #####


    ### --- LOG-LINEAR INSET GRAPH ---

    scatter!(average_saddle_index_densities_graph, 
        [-1.0 .+ (E0_bin_values_swap[non_zero_saddle_index_densities_indices_swap]./-solved_configuration_energy(cube)), 
        -1.0 .+ (E0_bin_values_slice[non_zero_saddle_index_densities_indices_slice]/-solved_configuration_energy(cube))],
        [log.(average_saddle_index_densities_swap[non_zero_saddle_index_densities_indices_swap]),
        log.(average_saddle_index_densities_slice_including_K_1_slice[non_zero_saddle_index_densities_indices_slice])]; 
        color=[alex_blue alex_orange], legend=false, inset=bbox(0.21,0.2,0.35,0.45), subplot=2,
        xlabel=L"E/|\!\!E_s|", ylabel=L"\log\langle k \rangle", yguidefontsize=12,xguidefontsize=12)
   

    ### --- SAVE AND DISPLAY GRAPH ---
    savefig(average_saddle_index_densities_graph, "results/final_paper_results/$(simulation_name)_average_saddle_index_densities.svg")
    display(average_saddle_index_densities_graph)
end